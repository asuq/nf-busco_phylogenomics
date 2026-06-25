#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." >/dev/null 2>&1 && pwd -P)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/nf-busco-frac-label.XXXXXX")"
OUT_DIR="$TMP_ROOT/output"
WORK_DIR="$TMP_ROOT/work"
NO_CONTAINER_CONFIG="$TMP_ROOT/no_container.config"
TRACE_FILE="$TMP_ROOT/trace.txt"

cleanup() {
    rm -rf "$TMP_ROOT"
}
trap cleanup EXIT

fail() {
    printf 'FAIL: %s\n' "$*" >&2
    exit 1
}

assert_path_exists() {
    [ -e "$1" ] || fail "expected path to exist: $1"
}

assert_trace_has_process() {
    grep -Fq "$1" "$TRACE_FILE" || fail "expected trace to include process: $1"
}

cd "$REPO_ROOT"

cat >"$NO_CONTAINER_CONFIG" <<'CONFIG'
params.max_cpus = 8
params.max_memory = 128.GB
params.max_time = 12.h
process.container = null
singularity.enabled = false
docker.enabled = false
apptainer.enabled = false
process {
    withLabel: process_single {
        cpus = 1
        memory = '1 GB'
        time = '1h'
    }
    withLabel: process_low {
        cpus = 1
        memory = '1 GB'
        time = '1h'
    }
    withLabel: process_medium {
        cpus = 1
        memory = '1 GB'
        time = '1h'
    }
    withLabel: process_high {
        cpus = 1
        memory = '1 GB'
        time = '1h'
    }
    withLabel: process_high_memory {
        cpus = 1
        memory = '1 GB'
        time = '1h'
    }
}
CONFIG

nextflow run main.nf \
    -c "$NO_CONTAINER_CONFIG" \
    -stub-run \
    -profile test,local \
    -ansi-log false \
    -work-dir "$WORK_DIR" \
    -with-trace "$TRACE_FILE" \
    --sample test/sample.csv \
    --lineage bacteria_odb10 \
    --outdir "$OUT_DIR"

assert_path_exists "$OUT_DIR/frac100pct_results/concat.faa"
assert_path_exists "$OUT_DIR/frac100pct_results/partitions.nex"
assert_path_exists "$OUT_DIR/frac100pct_results/frac100pct.iqtree"
assert_trace_has_process "concat_alignments"
assert_trace_has_process "run_iqtree"

printf 'Nextflow fraction publishDir stub test passed\n'
