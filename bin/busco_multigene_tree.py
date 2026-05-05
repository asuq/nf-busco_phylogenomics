#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    busco_multigene_tree.py <subcommand> [options]

Purpose:
    Create multi-gene phylogenomic tree from BUSCO single-copy orthologs.

Subcommands:
    all                         Run all steps: collect, select, align, infer
    collect                     Collect BUSCO FASTAs from directories or tar archives
    select                      Select shared genes & write gene lists
    align                       Align & trim gene alignments (mafft & trimAl)
    infer                       Concatenate & build tree (AMAS & IQ-TREE)

Options:
    -h, --help                  Show this
    -i, --input_dir             Path to the input directory which contains all BUSCO outputs
    -f, --fraction              Comma-spliced fractions for creating multi-gene phylogenetic tree
    -c, --cores                 Number of CPUs to execute [default: 8]
    -m, --mafft                 mafft options [default: --globalpair --maxiterate 1000 --thread $CORES]
    -t, --trimal                trimAl options
    -a, --amas                  AMAS options
    -q, --iqtree                IQ-TREE options
    -o, --out_dir               Output directory path [default: ./output]

Required Packages and Softwares:
    Biopython, tqdm, mafft, trimAl, AMAS, IQ-TREE

Author: Akito Shima (ASUQ)
Email: akito.shima@oist.jp
"""

import argparse
from decimal import Decimal, InvalidOperation, ROUND_CEILING
import io
import logging
import os
import shlex
import shutil
import subprocess
import sys
import tarfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import TextIO

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s")

SOFTWARES = ("mafft", "trimal", "AMAS.py", "iqtree")
BUSCO_SEQUENCE_DIRS = ("single_copy_busco_sequences", "multi_copy_busco_sequences")
BUSCO_SEQUENCE_ARCHIVES = tuple(f"{dirname}.tar.gz" for dirname in BUSCO_SEQUENCE_DIRS)

OPTION_DEFS = {
    "input_dir": dict(
        flags=("-i", "--input_dir"),
        kwargs=dict(type=Path, required=True, help="Path to BUSCO outputs"),
    ),
    "out_dir": dict(
        flags=("-o", "--out_dir"),
        kwargs=dict(
            type=Path,
            default=Path("./output"),
            help="Output directory [default: ./output]",
        ),
    ),
    "cores": dict(
        flags=("-c", "--cores"),
        kwargs=dict(type=int, default=8, help="CPUs to use [default: 8]"),
    ),
    "verbose": dict(
        flags=("-v", "--verbose"),
        kwargs=dict(action="store_true", help=argparse.SUPPRESS),
    ),
    "fraction": dict(
        flags=("-f", "--fraction"),
        kwargs=dict(
            type=str,
            default="0.9",
            help=(
                "Comma-spliced completeness fraction(s)\n"
                "(e.g. '0.8,0.99,0.999') [default: 0.9]"
            ),
        ),
    ),
    "mafft": dict(
        flags=("-m", "--mafft"),
        kwargs=dict(
            type=str,
            metavar="MAFFT_OPTIONS",
            default="--globalpair --maxiterate 1000 --thread $CORES",
            help=(
                "MAFFT options\n"
                "[default: mafft --globalpair --maxiterate 1000\n"
                "    --thread $CORES $INPUT > $OUTPUT]\n"
                "Parse with double quotes"
            ),
        ),
    ),
    "trimal": dict(
        flags=("-t", "--trimal"),
        kwargs=dict(
            type=str,
            metavar="TRIMAL_OPTIONS",
            default="-automated1",
            help="trimAl options\n"
            "[default: trimal -automated1 -in $INPUT -out $OUTPUT]\n"
            "Parse with double quotes",
        ),
    ),
    "amas": dict(
        flags=("-a", "--amas"),
        kwargs=dict(
            type=str,
            metavar="AMAS_OPTIONS",
            default="concat --in-format fasta --cores $CORES --data-type aa --part-format nexus",
            help="AMAS options\n"
            "[default: AMAS.py concat --in-files $INPUT --in-format fasta\n"
            "    --data-type aa --concat-out $OUT_concat.faa\n"
            "    --concat-part $OUT_partitions.txt\n"
            "    --part-format nexus --cores $Cores]\n"
            "Parse with double quotes",
        ),
    ),
    "iqtree": dict(
        flags=("-q", "--iqtree"),
        kwargs=dict(
            type=str,
            metavar="IQTREE_OPTIONS",
            default="-B 1000 -alrt 1000 -m MFP+MERGE -T AUTO",
            help="IQ-TREE options\n"
            "[default: iqtree -B 1000 -alrt 1000 -m MFP+MERGE\n"
            "    -T AUTO -s $INPUT -p $PARTITION -pre $PREFIX]\n"
            "Parse with double quotes",
        ),
    ),
}


SUBCMD_OPTS = {
    "collect": (
        ("input_dir", "out_dir", "cores", "verbose"),
        "Collect per-gene FASTA files from BUSCO outputs or BUSCO --tar archives",
    ),
    "select": (
        ("input_dir", "out_dir", "fraction", "cores", "verbose"),
        "Select shared genes & write gene lists",
    ),
    "align": (
        ("input_dir", "out_dir", "fraction", "cores", "mafft", "trimal", "verbose"),
        "Align & trim gene alignments (mafft & trimAl)",
    ),
    "infer": (
        ("input_dir", "out_dir", "fraction", "cores", "amas", "iqtree", "verbose"),
        "Concatenate & build tree (AMAS & IQ-TREE)",
    ),
    "all": (
        (
            "input_dir",
            "out_dir",
            "fraction",
            "cores",
            "mafft",
            "trimal",
            "amas",
            "iqtree",
            "verbose",
        ),
        "Run all steps: collect, select, align, infer",
    ),
}

# Python version check
REQUIRED = (3, 10)
if sys.version_info < REQUIRED:
    logging.fatal(
        f"this script requires Python {REQUIRED[0]}.{REQUIRED[1]}+, "
        f"but you’re running {sys.version_info.major}.{sys.version_info.minor}"
    )
    sys.exit(1)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build multi-gene phylogenomic tree from BUSCO single-copy orthologous gene",
        epilog="Required package and Softwares: Biopython, mafft, trimAl, AMAS, IQ-TREE",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    subs = parser.add_subparsers(dest="command", required=True)

    for cmd, (opt_keys, help) in SUBCMD_OPTS.items():
        p = subs.add_parser(
            cmd,
            help=help,
            formatter_class=argparse.RawTextHelpFormatter,
        )
        for key in opt_keys:
            flags = OPTION_DEFS[key]["flags"]
            kwargs = OPTION_DEFS[key]["kwargs"]
            p.add_argument(*flags, **kwargs)

    return parser.parse_args()


def check_software(tool: str) -> None:
    """Verify that a required software is available in PATH."""
    if shutil.which(tool) is None:
        logging.fatal(f"Required software '{tool}' not found in PATH.")
        sys.exit(1)


def run_cmd(cmd: list[str], stdout=None) -> None:
    """Execute a subprocess call, abort on failure"""
    try:
        logging.debug(f"Executing: {cmd}")
        subprocess.check_call(cmd, stdout=stdout)
    except subprocess.CalledProcessError as e:
        logging.fatal(f"Command failed {e.returncode}: {' '.join(cmd)}")
        sys.exit(e.returncode)


## Functions for parsing arguments
def format_decimal(value: Decimal) -> str:
    """Format a Decimal without trailing zeroes or scientific notation."""
    value_str = format(value.normalize(), "f")
    if "." in value_str:
        value_str = value_str.rstrip("0").rstrip(".")
    return value_str


def fraction_to_label(frac: Decimal) -> str:
    """Convert a fraction into an exact percentage label."""
    pct = frac * Decimal("100")
    return f"frac{format_decimal(pct).replace('.', 'p')}pct"


def parse_fractions(frac_str: str) -> list[Decimal]:
    """Parse comma-separated fractions into a sorted list of Decimal values."""
    fracs: set[Decimal] = set()
    for frac in frac_str.split(","):
        frac = frac.strip()
        try:
            frac = Decimal(frac)
        except InvalidOperation:
            logging.error("fractions must be comma-spliced numbers")
            sys.exit(1)
        if not 0 < frac <= 1:
            logging.error("fractions must be numbers between 0 and 1")
            sys.exit(1)

        fracs.add(frac.normalize())

    return sorted(fracs)


def load_genes_for_fraction(frac: Decimal, output_dir: Path) -> list[str]:
    """Load the list of selected genes for a given fraction from disk."""
    frac_label = fraction_to_label(frac)
    path = output_dir / f"{frac_label}_results/{frac_label}_genes.txt"
    with path.open() as f:
        lines = f.read().splitlines()
    return [line.strip() for line in lines[2:] if line.strip()]


def find_busco_sequence_sources(input_dir: Path) -> list[Path]:
    """Find BUSCO sequence directories and BUSCO --tar archives."""
    sources: list[Path] = []
    for root, dirs, files in os.walk(input_dir, followlinks=True):
        root_path = Path(root)
        for dirname in dirs:
            if dirname in BUSCO_SEQUENCE_DIRS:
                dirpath = root_path / dirname
                if dirpath.is_dir():
                    sources.append(dirpath)

        for filename in files:
            if filename in BUSCO_SEQUENCE_ARCHIVES:
                archive_path = root_path / filename
                if archive_path.is_file():
                    sources.append(archive_path)

    return sorted(sources)


def busco_sequence_source_name(source: Path) -> str:
    """Return the BUSCO sequence source name without a tar suffix."""
    if source.name.endswith(".tar.gz"):
        return source.name.removesuffix(".tar.gz")
    return source.name


def busco_sample_name(input_dir: Path, source: Path) -> str:
    """Extract the sample name from a BUSCO sequence source path."""
    try:
        return source.relative_to(input_dir).parts[-5]
    except (IndexError, ValueError):
        logging.error(f"Can't extract sample name from path: {source}")
        sys.exit(1)


def set_record_sample(record: SeqRecord, sample: str) -> SeqRecord:
    """Set BUSCO record identifiers to the sample name."""
    record.id = sample
    record.description = ""
    return record


def parse_busco_sequence_file(
    handle: TextIO,
    sample: str,
    first_only: bool,
) -> list[SeqRecord]:
    """Parse a BUSCO FASTA handle and return sample-labelled records."""
    if first_only:
        record = next(SeqIO.parse(handle, "fasta"), None)
        return [set_record_sample(record, sample)] if record else []

    return [set_record_sample(record, sample) for record in SeqIO.parse(handle, "fasta")]


def parse_busco_sequence_dir(
    seq_dir: Path,
    sample: str,
    first_only: bool,
) -> dict[str, list[SeqRecord]]:
    """Parse BUSCO FASTA files from an untarred sequence directory."""
    local: dict[str, list[SeqRecord]] = defaultdict(list)
    for faa_file in sorted(seq_dir.glob("*.faa")):
        logging.debug(f"Extracting gene seq from {str(faa_file)}")
        gene = faa_file.stem
        with faa_file.open() as handle:
            local[gene].extend(parse_busco_sequence_file(handle, sample, first_only))
    return local


def parse_busco_sequence_archive(
    archive_path: Path,
    sample: str,
    first_only: bool,
) -> dict[str, list[SeqRecord]]:
    """Stream BUSCO FASTA records from a BUSCO --tar archive."""
    local: dict[str, list[SeqRecord]] = defaultdict(list)
    with tarfile.open(archive_path, "r:gz") as archive:
        members = sorted(
            (
                member
                for member in archive.getmembers()
                if member.isfile() and Path(member.name).suffix == ".faa"
            ),
            key=lambda member: member.name,
        )
        for member in members:
            logging.debug(f"Extracting gene seq from {archive_path}:{member.name}")
            gene = Path(member.name).stem
            fasta_handle = archive.extractfile(member)
            if fasta_handle is None:
                continue
            with io.TextIOWrapper(fasta_handle, encoding="utf-8") as handle:
                local[gene].extend(parse_busco_sequence_file(handle, sample, first_only))
    return local


def parse_busco_sequence_source(
    source: Path,
    sample: str,
) -> dict[str, list[SeqRecord]]:
    """Parse BUSCO sequence records from a directory or tar archive."""
    source_name = busco_sequence_source_name(source)
    first_only = source_name == "multi_copy_busco_sequences"
    if source.is_dir():
        return parse_busco_sequence_dir(source, sample, first_only)
    return parse_busco_sequence_archive(source, sample, first_only)


def collect_gene_seqs(
    input_dir: Path,
    output_dir: Path,
    cores: int,
) -> tuple[dict[str, set[str]], set[str]]:
    """
    Parse single-copy BUSCO FASTAs with SeqIO, prefix record.id with sample name,
    and append into per-gene files.
    """

    raw_dir = output_dir / "seqs" / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    if any(raw_dir.iterdir()):
        logging.fatal(f"{raw_dir} is not empty — aborting to avoid mixing old results")
        sys.exit(1)

    logging.info("Collecting genes from BUSCO outputs")

    sequence_sources = find_busco_sequence_sources(input_dir)
    if not sequence_sources:
        logging.error(f"No BUSCO sequence sources found in {input_dir}")
        sys.exit(1)

    sample_sources: dict[str, list[Path]] = defaultdict(list)
    for source in sequence_sources:
        sample = busco_sample_name(input_dir, source)
        sample_sources[sample].append(source)

    logging.info(
        f"Found {len(sample_sources)} samples "
        f"(across {len(sequence_sources)} single/multi sources); parsing with {cores} threads"
    )

    def parse_sample(sample: str, sources: list[Path]):
        """Parse BUSCO sequence sources for a single sample."""
        local: dict[str, list[SeqRecord]] = defaultdict(list)
        for source in sources:
            for gene, records in parse_busco_sequence_source(source, sample).items():
                local[gene].extend(records)
        return sample, local

    results: list[tuple[str, dict[str, list[SeqRecord]]]] = []
    with ThreadPoolExecutor(max_workers=cores) as pool:
        futures = {
            pool.submit(parse_sample, sample, sources): sample
            for sample, sources in sample_sources.items()
        }
        for future in tqdm(as_completed(futures), total=len(futures), desc="Parsing samples"):
            try:
                results.append(future.result())
            except Exception as e:
                logging.error(f"Error parsing sample {futures[future]}: {e}")
                sys.exit(1)

    logging.info(f"Parsed {len(results)} samples successfully")

    gene_dict: dict[str, set[str]] = defaultdict(set)
    org_set: set[str] = set()
    gene_records: dict[str, list[SeqRecord]] = defaultdict(list)

    for org_name, local in results:
        if org_name in org_set:
            logging.error(f"Sample name is duplicated: {org_name}")
            sys.exit(1)
        org_set.add(org_name)

        for gene, recs in local.items():
            gene_dict[gene].add(org_name)
            gene_records[gene].extend(recs)

    def write_gene(item: tuple[str, list[SeqRecord]]) -> None:
        """Write out the gene records to a file"""
        gene, recs = item
        out_file = raw_dir / f"{gene}.faa"
        with out_file.open("w") as out:
            SeqIO.write(recs, out, "fasta")

    logging.info(f"Writing {len(gene_records)} gene files with {cores} threads")
    with ThreadPoolExecutor(max_workers=cores) as pool:
        list(
            tqdm(
                pool.map(write_gene, gene_records.items()),
                total=len(gene_records),
                desc="Writing gene files",
            )
        )

    logging.info(f"Finished writing {len(gene_records)} gene files")
    logging.debug(f"org_set: {org_set}")
    logging.debug(f"gene_dict: {gene_dict}")
    logging.info(f"Finished collecting gene sequences from {len(org_set)} samples")

    return gene_dict, org_set


def load_gene_dict(output_dir: Path) -> tuple[dict[str, set[str]], set[str]]:
    """Parse previously collected gene FASTAs to rebuild gene_dict and org_set."""
    gene_dict: dict[str, set[str]] = defaultdict(set)
    org_set: set[str] = set()

    raw_dir = output_dir / "seqs" / "raw"
    if not raw_dir.exists() or not any(raw_dir.iterdir()):
        logging.error(f"Raw directory '{raw_dir}' is missing or empty. Run 'collect' first.")
        sys.exit(1)

    logging.info("Loading gene sequences from raw directory")
    for faa_file in raw_dir.glob("*.faa"):
        gene = faa_file.stem
        for rec in SeqIO.parse(faa_file, "fasta"):
            org_set.add(rec.id)
            gene_dict[gene].add(rec.id)

    return gene_dict, org_set


def select_shared_genes(
    gene_dict: dict[str, set[str]], org_set: set[str], fractions: list[Decimal]
) -> dict[Decimal, list[str]]:
    """For each fraction, return list of genes present above the threshold."""

    logging.info("Selecting shared genes")

    total = len(org_set)
    frac_dict: dict[Decimal, list[str]] = {}
    for frac in fractions:
        threshold = int((Decimal(total) * frac).to_integral_value(rounding=ROUND_CEILING))
        frac_dict[frac] = [gene for gene, orgs in gene_dict.items() if len(orgs) >= threshold]

    logging.debug(f"frac_dict: {frac_dict}")
    return frac_dict


def write_gene_lists(frac_dict: dict[Decimal, list[str]], output_dir: Path) -> None:
    """Write out which genes pass completeness threshold."""

    logging.info("Writing out gene lists")

    for frac, genes in frac_dict.items():
        frac_label = fraction_to_label(frac)
        results_dir = output_dir / f"{frac_label}_results"
        results_dir.mkdir(parents=True, exist_ok=True)
        if any(results_dir.iterdir()):
            logging.fatal(f"{results_dir} is not empty — aborting to avoid mixing old results")
            sys.exit(1)

        file_path = results_dir / f"{frac_label}_genes.txt"
        with file_path.open("w") as out:
            out.write(f"Number of genes considered: {len(genes)}\n")
            out.write("Analyzed genes:\n")
            out.write("\n".join(genes) + "\n")


def align_and_trim(
    fractions: list[Decimal],
    output_dir: Path,
    mafft_opts: list[str],
    trimal_opts: list[str],
) -> None:
    """Align and trim all genes in the most inclusive set"""

    logging.info("Aligning the genes")

    seq_dir = output_dir / "seqs"
    raw_dir = seq_dir / "raw"
    aligned_dir = seq_dir / "aligned"
    trimmed_dir = seq_dir / "trimmed"

    for d in (aligned_dir, trimmed_dir):
        d.mkdir(parents=True, exist_ok=True)
        if any(d.iterdir()):
            logging.fatal(f"{d} is not empty — aborting to avoid mixing old results")
            sys.exit(1)

    smallest_frac = min(fractions)
    genes = load_genes_for_fraction(smallest_frac, output_dir)
    if not genes:
        logging.error(f"No genes for the smallest fraction {smallest_frac}.")
        sys.exit(1)

    logging.info(f"Aligning genes in fraction: {smallest_frac} which has {len(genes)} genes")

    for gene in tqdm(genes, desc="Align & trim"):
        infile = raw_dir / f"{gene}.faa"
        aligned = aligned_dir / f"{gene}_aligned.faa"
        trimmed = trimmed_dir / f"{gene}_trimmed.faa"

        logging.info(f"Running mafft {str(infile)} -> {str(aligned)}")
        with aligned.open("w") as out_f:
            run_cmd(["mafft"] + mafft_opts + [str(infile)], stdout=out_f)

        logging.info(f"Running trimAl {str(aligned)} -> {str(trimmed)}")
        run_cmd(["trimal"] + trimal_opts + ["-in", str(aligned), "-out", str(trimmed)])


def concat_alignments(
    frac_dict: dict[Decimal, list[str]], output_dir: Path, amas_opts: list[str]
) -> dict[Decimal, tuple[Path, Path]]:
    """Run AMAS to concatenate trimmed gene alignments using AMAS"""

    logging.info(f"Concatenating alignments")

    cafiles: dict[Decimal, tuple[Path, Path]] = {}
    for frac, genes in frac_dict.items():
        if not genes:
            logging.warning(f"No genes for fraction {frac}, skipping...")
            continue

        frac_label = fraction_to_label(frac)
        results_dir = output_dir / f"{frac_label}_results"
        concat_faa = results_dir / "concat.faa"
        partition_file = results_dir / "partitions.nex"

        trimmed_files = [str(output_dir / f"seqs/trimmed/{g}_trimmed.faa") for g in genes]

        logging.info(f"Running AMAS concat: fraction {frac}")
        cmd = (
            ["AMAS.py"]
            + amas_opts
            + [
                "--concat-out",
                str(concat_faa),
                "--concat-part",
                str(partition_file),
                "--in-files",
            ]
            + trimmed_files
        )
        run_cmd(cmd)
        cafiles[frac] = (concat_faa, partition_file)

    return cafiles


def run_iqtree(
    cafiles: dict[Decimal, tuple[Path, Path]], output_dir: Path, iqtree_opts: list[str]
) -> None:
    """Run IQ-TREE"""

    logging.info("Building trees")

    cafiles_sorted = {key: cafiles[key] for key in sorted(cafiles.keys(), reverse=True)}

    for frac, (concat_faa, partition_file) in cafiles_sorted.items():
        frac_label = fraction_to_label(frac)
        prefix = output_dir / f"{frac_label}_results/{frac_label}"

        logging.info(f"Running IQ-TREE on {str(concat_faa)} and {str(partition_file)}")
        run_cmd(
            ["iqtree"]
            + iqtree_opts
            + ["-s", str(concat_faa), "-p", str(partition_file), "-pre", str(prefix)]
        )


def main() -> None:
    # Verify software in the PATH
    for tool in SOFTWARES:
        check_software(tool)

    args = parse_args()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logging.debug(f"parsed arguments: {args!r}")

    if not args.input_dir.is_dir():
        logging.error(f"Input {str(args.input_dir)} is not directory")
        sys.exit(1)

    logging.info("Starting process...")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.command == "collect":
        logging.info("Running subcommand: collect")
        gene_dict, org_set = collect_gene_seqs(args.input_dir, args.out_dir, args.cores)

    elif args.command == "select":
        logging.info("Running subcommand: select")
        fracs = parse_fractions(args.fraction)
        gene_dict, org_set = load_gene_dict(args.out_dir)
        frac_dict = select_shared_genes(gene_dict, org_set, fracs)
        write_gene_lists(frac_dict, args.out_dir)

    elif args.command == "align":
        logging.info("Running subcommand: align")
        fracs = parse_fractions(args.fraction)
        args.mafft = shlex.split(args.mafft.replace("$CORES", str(args.cores)))
        args.trimal = shlex.split(args.trimal)

        align_and_trim(fracs, args.out_dir, args.mafft, args.trimal)

    elif args.command == "infer":
        logging.info("Running subcommand: infer")
        fracs = parse_fractions(args.fraction)
        args.amas = shlex.split(args.amas.replace("$CORES", str(args.cores)))
        args.iqtree = shlex.split(args.iqtree)

        frac_dict = {}
        for frac in fracs:
            genes = load_genes_for_fraction(frac, args.out_dir)
            frac_dict[frac] = genes
        cafiles = concat_alignments(frac_dict, args.out_dir, args.amas)
        run_iqtree(cafiles, args.out_dir, args.iqtree)

    elif args.command == "all":
        logging.info("Running subcommand: all")
        fracs = parse_fractions(args.fraction)
        args.mafft = shlex.split(args.mafft.replace("$CORES", str(args.cores)))
        args.trimal = shlex.split(args.trimal)
        args.amas = shlex.split(args.amas.replace("$CORES", str(args.cores)))
        args.iqtree = shlex.split(args.iqtree)

        gene_dict, org_set = collect_gene_seqs(args.input_dir, args.out_dir, args.cores)
        frac_dict = select_shared_genes(gene_dict, org_set, fracs)
        write_gene_lists(frac_dict, args.out_dir)
        align_and_trim(fracs, args.out_dir, args.mafft, args.trimal)
        cafiles = concat_alignments(frac_dict, args.out_dir, args.amas)
        run_iqtree(cafiles, args.out_dir, args.iqtree)

    else:
        logging.error(f"Unknown command {args.command}")
        sys.exit(1)

    logging.info("Job completed")


if __name__ == "__main__":
    main()
