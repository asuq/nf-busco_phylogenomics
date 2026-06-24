"""Tests for building representative BUSCO phylogenomics samplesheets."""

from __future__ import annotations

import csv
import importlib.util
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "helpers"
    / "build_species_representative_sample.py"
)
SPEC = importlib.util.spec_from_file_location(
    "build_species_representative_sample",
    MODULE_PATH,
)
builder = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = builder
SPEC.loader.exec_module(builder)


def write_text(path: Path, text: str) -> None:
    """Write text after creating parent directories."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV file into dictionaries."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


class BuildSpeciesRepresentativeSampleTests(unittest.TestCase):
    """Exercise representative samplesheet conversion behaviour."""

    def assert_builder_fails(self, argv: list[str]) -> None:
        """Assert that the helper exits with an error for invalid inputs."""
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit) as context:
                builder.main(argv)
        self.assertEqual(context.exception.code, 1)

    def test_builds_sorted_representative_samplesheet(self) -> None:
        """Representatives are selected, sanitised, and sorted by cluster/accession."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            master_tsv = root / "master.tsv"
            sample_csv = root / "sample.csv"
            output_csv = root / "species_representatives.csv"
            write_text(root / "genomes" / "gca_002.fa", ">b\nAC\n")
            write_text(root / "genomes" / "gca_001.fa", ">a\nAC\n")
            write_text(root / "genomes" / "gca_003.fa", ">c\nAC\n")
            write_text(
                master_tsv,
                (
                    "Cluster_ID\taccession\tOrganism_Name\tIs_Representative\n"
                    "cluster-2\tGCA_002.1\tBeta species\tno\n"
                    "cluster-2\tGCA_002.2\tBeta species representative\tyes\n"
                    "cluster-1\tGCA_001.1\tAlpha/species\tyes\n"
                ),
            )
            write_text(
                sample_csv,
                (
                    "accession,is_new,assembly_level,genome_fasta\n"
                    "GCA_002.2,false,Complete Genome,genomes/gca_002.fa\n"
                    "GCA_001.1,true,Scaffold,genomes/gca_001.fa\n"
                    "GCA_003.1,true,Scaffold,genomes/gca_003.fa\n"
                ),
            )

            self.assertEqual(
                builder.main([str(master_tsv), str(sample_csv), str(output_csv)]),
                0,
            )

            self.assertEqual(
                read_csv(output_csv),
                [
                    {
                        "sample": "cluster_1_GCA_001_1_Alpha_species",
                        "fasta": "genomes/gca_001.fa",
                    },
                    {
                        "sample": "cluster_2_GCA_002_2_Beta_species_representative",
                        "fasta": "genomes/gca_002.fa",
                    },
                ],
            )

    def test_writes_absolute_paths_when_requested(self) -> None:
        """Absolute mode resolves relative FASTA paths from the sample CSV directory."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fasta = root / "genomes" / "gca_001.fa"
            master_tsv = root / "master.tsv"
            sample_csv = root / "sample.csv"
            output_csv = root / "output.csv"
            write_text(fasta, ">a\nAC\n")
            write_text(
                master_tsv,
                (
                    "Cluster_ID\taccession\tOrganism_Name\tIs_Representative\n"
                    "c1\tGCA_001.1\tAlpha\tyes\n"
                ),
            )
            write_text(
                sample_csv,
                "accession,genome_fasta\nGCA_001.1,genomes/gca_001.fa\n",
            )

            self.assertEqual(
                builder.main(
                    [
                        str(master_tsv),
                        str(sample_csv),
                        str(output_csv),
                        "--path-mode",
                        "absolute",
                    ]
                ),
                0,
            )

            self.assertEqual(read_csv(output_csv)[0]["fasta"], str(fasta.resolve()))

    def test_accepts_configurable_columns(self) -> None:
        """Column override flags support non-default nf-annotation exports."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            master_tsv = root / "master.tsv"
            sample_csv = root / "sample.csv"
            output_csv = root / "output.csv"
            write_text(root / "genome.fa", ">a\nAC\n")
            write_text(
                master_tsv,
                "cluster\tassembly\tname\tchosen\nA-1\tASM1\tCustom name\ttrue\n",
            )
            write_text(sample_csv, "assembly,fasta\nASM1,genome.fa\n")

            self.assertEqual(
                builder.main(
                    [
                        str(master_tsv),
                        str(sample_csv),
                        str(output_csv),
                        "--cluster-column",
                        "cluster",
                        "--accession-column",
                        "assembly",
                        "--organism-column",
                        "name",
                        "--representative-column",
                        "chosen",
                        "--representative-value",
                        "true",
                        "--sample-accession-column",
                        "assembly",
                        "--fasta-column",
                        "fasta",
                    ]
                ),
                0,
            )

            self.assertEqual(
                read_csv(output_csv),
                [{"sample": "A_1_ASM1_Custom_name", "fasta": "genome.fa"}],
            )

    def test_fails_on_duplicate_representative_accession(self) -> None:
        """Duplicate representative accessions in the master TSV fail."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            master_tsv = root / "master.tsv"
            sample_csv = root / "sample.csv"
            output_csv = root / "output.csv"
            write_text(root / "genome.fa", ">a\nAC\n")
            write_text(
                master_tsv,
                (
                    "Cluster_ID\taccession\tOrganism_Name\tIs_Representative\n"
                    "c1\tGCA_001.1\tAlpha\tyes\n"
                    "c2\tGCA_001.1\tAlpha duplicate\tyes\n"
                ),
            )
            write_text(sample_csv, "accession,genome_fasta\nGCA_001.1,genome.fa\n")

            self.assert_builder_fails(
                [str(master_tsv), str(sample_csv), str(output_csv)]
            )

    def test_fails_when_representative_is_missing_from_sample_csv(self) -> None:
        """Representative accessions absent from the sample CSV fail."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            master_tsv = root / "master.tsv"
            sample_csv = root / "sample.csv"
            output_csv = root / "output.csv"
            write_text(root / "other.fa", ">b\nAC\n")
            write_text(
                master_tsv,
                (
                    "Cluster_ID\taccession\tOrganism_Name\tIs_Representative\n"
                    "c1\tGCA_001.1\tAlpha\tyes\n"
                ),
            )
            write_text(sample_csv, "accession,genome_fasta\nGCA_002.1,other.fa\n")

            self.assert_builder_fails(
                [str(master_tsv), str(sample_csv), str(output_csv)]
            )

    def test_fails_when_fasta_path_is_missing(self) -> None:
        """Missing FASTA files fail by default."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            master_tsv = root / "master.tsv"
            sample_csv = root / "sample.csv"
            output_csv = root / "output.csv"
            write_text(
                master_tsv,
                (
                    "Cluster_ID\taccession\tOrganism_Name\tIs_Representative\n"
                    "c1\tGCA_001.1\tAlpha\tyes\n"
                ),
            )
            write_text(sample_csv, "accession,genome_fasta\nGCA_001.1,missing.fa\n")

            self.assert_builder_fails(
                [str(master_tsv), str(sample_csv), str(output_csv)]
            )

    def test_skip_fasta_check_allows_missing_paths(self) -> None:
        """Missing FASTA files are allowed when checking is explicitly skipped."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            master_tsv = root / "master.tsv"
            sample_csv = root / "sample.csv"
            output_csv = root / "output.csv"
            write_text(
                master_tsv,
                (
                    "Cluster_ID\taccession\tOrganism_Name\tIs_Representative\n"
                    "c1\tGCA_001.1\tAlpha\tyes\n"
                ),
            )
            write_text(sample_csv, "accession,genome_fasta\nGCA_001.1,missing.fa\n")

            self.assertEqual(
                builder.main(
                    [
                        str(master_tsv),
                        str(sample_csv),
                        str(output_csv),
                        "--skip-fasta-check",
                    ]
                ),
                0,
            )

            self.assertEqual(read_csv(output_csv)[0]["fasta"], "missing.fa")


if __name__ == "__main__":
    unittest.main()
