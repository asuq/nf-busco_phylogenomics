"""Tests for collecting BUSCO genes from tarred and untarred outputs."""

import importlib.util
import io
import tarfile
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


MODULE_PATH = Path(__file__).resolve().parents[1] / "bin" / "busco_multigene_tree.py"
SPEC = importlib.util.spec_from_file_location("busco_multigene_tree", MODULE_PATH)
busco_multigene_tree = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(busco_multigene_tree)


def busco_sequences_dir(input_dir: Path, sample: str) -> Path:
    """Return a BUSCO-like sequence output directory for a sample."""
    return input_dir / sample / "busco_output" / "run_test_odb10" / "busco_sequences"


def write_text(path: Path, text: str) -> None:
    """Write text to a file after creating its parent directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def write_tar(path: Path, members: dict[str, str]) -> None:
    """Write text members to a gzipped tar archive."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(path, "w:gz") as archive:
        for name, text in members.items():
            data = text.encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))


def read_records(path: Path) -> list[SeqRecord]:
    """Read FASTA records from a path."""
    with path.open() as handle:
        return list(SeqIO.parse(handle, "fasta"))


class BuscoTarCollectionTests(unittest.TestCase):
    """Exercise BUSCO collection across tarred and untarred layouts."""

    def test_collects_untarred_busco_directories(self) -> None:
        """Directory BUSCO outputs still collect single and multi-copy genes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_dir = root / "input"
            output_dir = root / "output"
            seq_dir = busco_sequences_dir(input_dir, "sample_a")

            write_text(
                seq_dir / "single_copy_busco_sequences" / "geneA.faa",
                ">geneA_a\nMAAA\n",
            )
            write_text(
                seq_dir / "multi_copy_busco_sequences" / "geneB.faa",
                ">geneB_first\nMBBB\n>geneB_second\nMCCC\n",
            )

            gene_dict, org_set = busco_multigene_tree.collect_gene_seqs(
                input_dir, output_dir, 1
            )

            self.assertEqual(org_set, {"sample_a"})
            self.assertEqual(gene_dict["geneA"], {"sample_a"})
            self.assertEqual(gene_dict["geneB"], {"sample_a"})
            gene_b_records = read_records(output_dir / "seqs" / "raw" / "geneB.faa")
            self.assertEqual(len(gene_b_records), 1)
            self.assertEqual(gene_b_records[0].id, "sample_a")
            self.assertEqual(str(gene_b_records[0].seq), "MBBB")

    def test_collects_tarred_busco_archives(self) -> None:
        """BUSCO --tar archives collect without extraction to the output tree."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_dir = root / "input"
            output_dir = root / "output"
            seq_dir = busco_sequences_dir(input_dir, "sample_a")

            write_tar(
                seq_dir / "single_copy_busco_sequences.tar.gz",
                {"single_copy_busco_sequences/geneA.faa": ">geneA_a\nMAAA\n"},
            )
            write_tar(
                seq_dir / "multi_copy_busco_sequences.tar.gz",
                {
                    "multi_copy_busco_sequences/geneB.faa": (
                        ">geneB_first\nMBBB\n>geneB_second\nMCCC\n"
                    )
                },
            )

            gene_dict, org_set = busco_multigene_tree.collect_gene_seqs(
                input_dir, output_dir, 1
            )

            self.assertEqual(org_set, {"sample_a"})
            self.assertEqual(gene_dict["geneA"], {"sample_a"})
            self.assertEqual(gene_dict["geneB"], {"sample_a"})
            gene_a_records = read_records(output_dir / "seqs" / "raw" / "geneA.faa")
            self.assertEqual(gene_a_records[0].id, "sample_a")
            self.assertEqual(str(gene_a_records[0].seq), "MAAA")

    def test_collects_mixed_tarred_and_untarred_outputs(self) -> None:
        """Mixed BUSCO layouts across samples collect into the same gene file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_dir = root / "input"
            output_dir = root / "output"

            sample_a_seq_dir = busco_sequences_dir(input_dir, "sample_a")
            write_text(
                sample_a_seq_dir / "single_copy_busco_sequences" / "geneA.faa",
                ">geneA_a\nMAAA\n",
            )

            sample_b_seq_dir = busco_sequences_dir(input_dir, "sample_b")
            write_tar(
                sample_b_seq_dir / "single_copy_busco_sequences.tar.gz",
                {"single_copy_busco_sequences/geneA.faa": ">geneA_b\nMBBB\n"},
            )

            gene_dict, org_set = busco_multigene_tree.collect_gene_seqs(
                input_dir, output_dir, 1
            )

            self.assertEqual(org_set, {"sample_a", "sample_b"})
            self.assertEqual(gene_dict["geneA"], {"sample_a", "sample_b"})
            records = read_records(output_dir / "seqs" / "raw" / "geneA.faa")
            self.assertEqual({record.id for record in records}, {"sample_a", "sample_b"})
            self.assertEqual(
                {record.id: str(record.seq) for record in records},
                {"sample_a": "MAAA", "sample_b": "MBBB"},
            )

    def test_multi_copy_archive_uses_first_record(self) -> None:
        """Tarred multi-copy BUSCO genes keep only the first FASTA record."""
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            input_dir = root / "input"
            output_dir = root / "output"
            seq_dir = busco_sequences_dir(input_dir, "sample_a")

            write_tar(
                seq_dir / "multi_copy_busco_sequences.tar.gz",
                {
                    "multi_copy_busco_sequences/geneB.faa": (
                        ">geneB_first\nMBBB\n>geneB_second\nMCCC\n"
                    )
                },
            )

            busco_multigene_tree.collect_gene_seqs(input_dir, output_dir, 1)

            records = read_records(output_dir / "seqs" / "raw" / "geneB.faa")
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].id, "sample_a")
            self.assertEqual(str(records[0].seq), "MBBB")


if __name__ == "__main__":
    unittest.main()
