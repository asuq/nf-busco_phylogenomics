#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Usage:
    busco_multigene_tree.py <subcommand> [options]

Purpose:
    Create multi-gene phylogenomic tree from BUSCO single-copy orthologs.

Subcommands:
    all                         Run all steps: collect, select, align, infer
    collect                     Collect per-gene FASTA files from BUSCO outputs
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
import fnmatch
import logging
import math
import os
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from concurrent.futures import Future, ThreadPoolExecutor, as_completed
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm


logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s: %(message)s")

SOFTWARES = ("mafft", "trimal", "AMAS.py", "iqtree")

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
            help=f"Comma-spliced completeness fraction(s)\n(e.g. '0.8,0.9') [default: 0.9]",
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
        "Collect per-gene FASTA files from BUSCO outputs",
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
def parse_fractions(frac_str: str) -> list[float]:
    """Parse comma-separated fractions into a sorted list of floats."""
    fracs = set()
    for frac in frac_str.split(","):
        try:
            frac = float(frac)
        except ValueError:
            logging.error("fractions must be comma-spliced numbers")
            sys.exit(1)
        if not 0 < frac <= 1:
            logging.error("fractions must be numbers between 0 and 1")
            sys.exit(1)

        fracs.add(round(frac, 2))

    return sorted(fracs)


def load_genes_for_fraction(frac: float, output_dir: Path) -> list[str]:
    """Load the list of selected genes for a given fraction from disk."""
    pct = int(frac * 100)
    path = output_dir / f"frac{pct}pct_results/frac{pct}pct_genes.txt"
    with path.open() as f:
        lines = f.read().splitlines()
    return [line.strip() for line in lines[2:] if line.strip()]


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

    # 1. find both single- and multi-copy directories
    # Note: recurse_symlinks requires Python 3.13+
    # all_dirs: list[Path] = [
    #     d
    #     for d in input_dir.rglob("*_copy_busco_sequences", recurse_symlinks=True)
    #     if d.is_dir()
    # ]

    all_dirs: list[Path] = []
    for root, dirs, _ in os.walk(input_dir, followlinks=True):
        for dirname in dirs:
            if fnmatch.fnmatch(dirname, "*_copy_busco_sequences"):
                dirpath = Path(root) / dirname
                if dirpath.is_dir():
                    all_dirs.append(dirpath)

    if not all_dirs:
        logging.error(f"No BUSCO output directories found in {input_dir}")
        sys.exit(1)

    # 2. group by sample name (extracted the same way you did before)
    sample_dirs: dict[str, list[Path]] = defaultdict(list)
    for d in all_dirs:
        try:
            sample = d.relative_to(input_dir).parts[-5]
        except IndexError:
            logging.error(f"Can't extract sample name from path: {d}")
            sys.exit(1)
        sample_dirs[sample].append(d)

    logging.info(
        f"Found {len(sample_dirs)} samples "
        f"(across {len(all_dirs)} single/multi dirs); parsing with {cores} threads"
    )

    def parse_sample(sample: str, dirs: list[Path]):
        """Function to parse a single sample directory"""

        local: dict[str, list[SeqRecord]] = defaultdict(list)
        for seq_dir in dirs:
            is_multi = seq_dir.name == "multi_copy_busco_sequences"
            for faa_file in seq_dir.glob("*.faa"):
                logging.debug(f"Extracting gene seq from {str(faa_file)}")
                gene = faa_file.stem
                if is_multi:
                    # only first record
                    rec = next(SeqIO.parse(faa_file, "fasta"), None)
                    if rec:
                        rec.id = sample
                        rec.description = ""
                        local[gene].append(rec)
                else:
                    # all (should be one) record(s)
                    for rec in SeqIO.parse(faa_file, "fasta"):
                        rec.id = sample
                        rec.description = ""
                        local[gene].append(rec)
        return sample, local

    # 4. run in parallel over unique samples
    results: list[tuple[str, dict[str, list[SeqRecord]]]] = []
    with ThreadPoolExecutor(max_workers=cores) as pool:
        futures = {
            pool.submit(parse_sample, sample, dirs): sample
            for sample, dirs in sample_dirs.items()
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
    gene_dict: dict[str, set[str]], org_set: set[str], fractions: list[float]
) -> dict[float, list[str]]:
    """For each fraction, return list of genes present above the threshold."""

    logging.info("Selecting shared genes")

    total = len(org_set)
    frac_dict: dict[float, list[str]] = {}
    for frac in fractions:
        threshold = math.ceil(total * frac)
        frac_dict[frac] = [gene for gene, orgs in gene_dict.items() if len(orgs) >= threshold]

    logging.debug(f"frac_dict: {frac_dict}")
    return frac_dict


def write_gene_lists(frac_dict: dict[float, list[str]], output_dir: Path) -> None:
    """Write out which genes pass completeness threshold."""

    logging.info("Writing out gene lists")

    for frac, genes in frac_dict.items():
        pct = int(frac * 100)
        results_dir = output_dir / f"frac{pct}pct_results"
        results_dir.mkdir(parents=True, exist_ok=True)
        if any(results_dir.iterdir()):
            logging.fatal(f"{results_dir} is not empty — aborting to avoid mixing old results")
            sys.exit(1)

        file_path = results_dir / f"frac{pct}pct_genes.txt"
        with file_path.open("w") as out:
            out.write(f"Number of genes considered: {len(genes)}\n")
            out.write("Analyzed genes:\n")
            out.write("\n".join(genes) + "\n")


def align_and_trim(
    fractions: list[float],
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
    frac_dict: dict[float, list[str]], output_dir: Path, amas_opts: list[str]
) -> dict[float, tuple[Path, Path]]:
    """Run AMAS to concatenate trimmed gene alignments using AMAS"""

    logging.info(f"Concatenating alignments")

    cafiles: dict[float, tuple[Path, Path]] = {}
    for frac, genes in frac_dict.items():
        if not genes:
            logging.warning(f"No genes for fraction {frac}, skipping...")
            continue

        pct = int(frac * 100)
        results_dir = output_dir / f"frac{pct}pct_results"
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
    cafiles: dict[float, tuple[Path, Path]], output_dir: Path, iqtree_opts: list[str]
) -> None:
    """Run IQ-TREE"""

    logging.info("Building trees")

    cafiles_sorted = {key: cafiles[key] for key in sorted(cafiles.keys(), reverse=True)}

    for frac, (concat_faa, partition_file) in cafiles_sorted.items():
        pct = int(frac * 100)

        prefix = output_dir / f"frac{pct}pct_results/frac{pct}pct"

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
