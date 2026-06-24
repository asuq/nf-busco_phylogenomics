#!/usr/bin/env python3
"""Build a BUSCO phylogenomics samplesheet from nf-annotation outputs."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class ColumnConfig:
    """Store input column names and representative filter values."""

    cluster: str
    accession: str
    organism: str
    representative: str
    representative_value: str
    sample_accession: str
    fasta: str


@dataclass(frozen=True)
class Representative:
    """Store representative metadata needed to build a sample name."""

    cluster: str
    accession: str
    organism: str


@dataclass(frozen=True)
class SampleRow:
    """Store one output samplesheet row and its stable sort keys."""

    cluster: str
    accession: str
    sample: str
    fasta: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Build a BUSCO phylogenomics sample CSV from an nf-annotation "
            "master TSV and nf-annotation sample CSV."
        )
    )
    parser.add_argument("master_tsv", type=Path, help="nf-annotation master TSV")
    parser.add_argument(
        "sample_csv",
        type=Path,
        help="nf-annotation sample CSV containing FASTA paths",
    )
    parser.add_argument("output_csv", type=Path, help="Output sample CSV")
    parser.add_argument(
        "--cluster-column",
        default="Cluster_ID",
        help="Master TSV cluster column [default: Cluster_ID]",
    )
    parser.add_argument(
        "--accession-column",
        default="accession",
        help="Master TSV accession column [default: accession]",
    )
    parser.add_argument(
        "--organism-column",
        default="Organism_Name",
        help="Master TSV organism column [default: Organism_Name]",
    )
    parser.add_argument(
        "--representative-column",
        default="Is_Representative",
        help="Master TSV representative flag column [default: Is_Representative]",
    )
    parser.add_argument(
        "--representative-value",
        default="yes",
        help="Value marking representative rows [default: yes]",
    )
    parser.add_argument(
        "--sample-accession-column",
        default="accession",
        help="Sample CSV accession column [default: accession]",
    )
    parser.add_argument(
        "--fasta-column",
        default="genome_fasta",
        help="Sample CSV FASTA path column [default: genome_fasta]",
    )
    parser.add_argument(
        "--path-mode",
        choices=("preserve", "absolute"),
        default="preserve",
        help=(
            "Write FASTA paths as preserved input text or absolute paths "
            "[default: preserve]"
        ),
    )
    parser.add_argument(
        "--skip-fasta-check",
        action="store_true",
        help="Do not check that FASTA paths exist",
    )
    return parser.parse_args(argv)


def fail(message: str) -> None:
    """Exit with a clear error message."""
    print(f"ERROR: {message}", file=sys.stderr)
    raise SystemExit(1)


def require_columns(
    fieldnames: list[str] | None,
    required: Sequence[str],
    label: str,
) -> None:
    """Require that a delimited input contains named columns."""
    if fieldnames is None:
        fail(f"{label} has no header")

    missing = [name for name in required if name not in fieldnames]
    if missing:
        fail(f"{label} missing required column(s): {', '.join(missing)}")


def sanitise_sample_part(value: str, label: str) -> str:
    """Return a filesystem-friendly sample-name component."""
    stripped = value.strip()
    if not stripped:
        fail(f"representative row has empty {label}")

    sample_part = re.sub(r"[^A-Za-z0-9]+", "_", stripped)
    sample_part = re.sub(r"^_+|_+$", "", sample_part)
    if not sample_part:
        fail(f"representative row has unusable {label}: {value}")

    return sample_part


def read_representatives(
    master_tsv: Path,
    columns: ColumnConfig,
) -> dict[str, Representative]:
    """Read representative accessions and sample-name metadata from the master TSV."""
    required = (
        columns.cluster,
        columns.accession,
        columns.organism,
        columns.representative,
    )
    representatives: dict[str, Representative] = {}

    with master_tsv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames, required, f"TSV {master_tsv}")

        for row in reader:
            representative_value = (row[columns.representative] or "").strip()
            if representative_value != columns.representative_value:
                continue

            accession = (row[columns.accession] or "").strip()
            if not accession:
                fail("found representative row with empty accession in master TSV")
            if accession in representatives:
                fail(f"duplicate representative accession in master TSV: {accession}")

            representatives[accession] = Representative(
                cluster=sanitise_sample_part(
                    row[columns.cluster] or "",
                    columns.cluster,
                ),
                accession=sanitise_sample_part(accession, columns.accession),
                organism=sanitise_sample_part(
                    row[columns.organism] or "",
                    columns.organism,
                ),
            )

    if not representatives:
        fail("no representative rows found in master TSV")

    return representatives


def existing_fasta_path(fasta_text: str, sample_csv: Path) -> Path | None:
    """Return the first existing FASTA path candidate for a sample CSV value."""
    raw_path = Path(fasta_text).expanduser()
    candidates = [raw_path]
    if not raw_path.is_absolute():
        candidates.append(sample_csv.parent / raw_path)

    for candidate in candidates:
        if candidate.is_file():
            return candidate

    return None


def output_fasta_path(
    fasta_text: str,
    sample_csv: Path,
    accession: str,
    path_mode: str,
    skip_fasta_check: bool,
) -> str:
    """Validate and format a FASTA path for the output samplesheet."""
    if not fasta_text.strip():
        fail(f"empty FASTA path for accession: {accession}")

    if skip_fasta_check:
        if path_mode == "absolute":
            raw_path = Path(fasta_text).expanduser()
            if not raw_path.is_absolute():
                raw_path = sample_csv.parent / raw_path
            return str(raw_path.resolve(strict=False))
        return fasta_text

    fasta_path = existing_fasta_path(fasta_text, sample_csv)
    if fasta_path is None:
        fail(f"FASTA file not found for accession {accession}: {fasta_text}")

    if path_mode == "absolute":
        return str(fasta_path.resolve())

    return fasta_text


def read_fasta_paths(
    sample_csv: Path,
    wanted_accessions: set[str],
    columns: ColumnConfig,
    path_mode: str,
    skip_fasta_check: bool,
) -> dict[str, str]:
    """Read FASTA paths for representative accessions from the sample CSV."""
    required = (columns.sample_accession, columns.fasta)
    fasta_paths: dict[str, str] = {}
    seen_accessions: set[str] = set()

    with sample_csv.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        require_columns(reader.fieldnames, required, f"CSV {sample_csv}")

        for row in reader:
            accession = (row[columns.sample_accession] or "").strip()
            if accession not in wanted_accessions:
                continue

            if accession in seen_accessions:
                fail(f"duplicate accession in sample CSV: {accession}")
            seen_accessions.add(accession)

            fasta_paths[accession] = output_fasta_path(
                fasta_text=row[columns.fasta] or "",
                sample_csv=sample_csv,
                accession=accession,
                path_mode=path_mode,
                skip_fasta_check=skip_fasta_check,
            )

    missing = sorted(wanted_accessions - fasta_paths.keys())
    if missing:
        preview = ", ".join(missing[:20])
        suffix = "" if len(missing) <= 20 else f" ... and {len(missing) - 20} more"
        fail(f"representative accession(s) missing from sample CSV: {preview}{suffix}")

    return fasta_paths


def build_rows(
    representatives: dict[str, Representative],
    fasta_paths: dict[str, str],
) -> list[SampleRow]:
    """Build sorted samplesheet rows from representative metadata and FASTA paths."""
    rows = [
        SampleRow(
            cluster=representative.cluster,
            accession=representative.accession,
            sample=(
                f"{representative.cluster}_"
                f"{representative.accession}_"
                f"{representative.organism}"
            ),
            fasta=fasta_paths[raw_accession],
        )
        for raw_accession, representative in representatives.items()
    ]
    return sorted(rows, key=lambda row: (row.cluster, row.accession))


def write_output(output_csv: Path, rows: Sequence[SampleRow]) -> None:
    """Write BUSCO phylogenomics sample CSV rows."""
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sample", "fasta"])
        for row in rows:
            writer.writerow([row.sample, row.fasta])


def build_column_config(args: argparse.Namespace) -> ColumnConfig:
    """Build a column configuration from parsed arguments."""
    return ColumnConfig(
        cluster=args.cluster_column,
        accession=args.accession_column,
        organism=args.organism_column,
        representative=args.representative_column,
        representative_value=args.representative_value,
        sample_accession=args.sample_accession_column,
        fasta=args.fasta_column,
    )


def main(argv: Sequence[str] | None = None) -> int:
    """Run the representative samplesheet builder."""
    args = parse_args(argv)

    if not args.master_tsv.is_file():
        fail(f"master TSV not found: {args.master_tsv}")
    if not args.sample_csv.is_file():
        fail(f"sample CSV not found: {args.sample_csv}")

    columns = build_column_config(args)
    representatives = read_representatives(args.master_tsv, columns)
    fasta_paths = read_fasta_paths(
        sample_csv=args.sample_csv,
        wanted_accessions=set(representatives),
        columns=columns,
        path_mode=args.path_mode,
        skip_fasta_check=args.skip_fasta_check,
    )
    rows = build_rows(representatives, fasta_paths)
    write_output(args.output_csv, rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
