#!/usr/bin/env bash
# End-to-end smoke test on a small public dataset.
#
# Downloads the 10x Genomics PBMC 3k dataset (~8 MB, 2,700 peripheral blood mononuclear
# cells) and runs the full pipeline: QC -> doublets -> clustering -> markers -> report.
#
# A second sample is created by subsampling 50% of the cells and writing it as .h5ad.
# It is NOT a biological replicate - it exists so the test exercises the multi-sample
# paths (concatenation, per-sample doublet scoring, cluster composition, Harmony
# integration) and both the 10x-directory and .h5ad readers in one run.
#
# Usage:
#   ./test/run_test.sh              # download if needed, then run
#   ./test/run_test.sh --dry-run    # build the DAG only, run nothing
#
# Runtime is roughly 1-2 minutes on a laptop.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
REPO_ROOT="$PWD"

DRY_RUN=""
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN="-n"

# Snakemake manages the conda envs by default. Set USE_CONDA=0 if scanpy and R are
# already on PATH.
CONDA_FLAG="--use-conda"
[[ "${USE_CONDA:-1}" == "0" ]] && CONDA_FLAG=""

DATA_DIR="test/data"
mkdir -p "$DATA_DIR"

URL="https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
MTX_DIR="$DATA_DIR/filtered_gene_bc_matrices/hg19"

if [[ ! -d "$MTX_DIR" ]]; then
    echo "==> Downloading PBMC 3k (~8 MB)"
    curl -sfL -o "$DATA_DIR/pbmc3k.tar.gz" "$URL"
    tar xzf "$DATA_DIR/pbmc3k.tar.gz" -C "$DATA_DIR"
    rm -f "$DATA_DIR/pbmc3k.tar.gz"
else
    echo "==> PBMC 3k already present"
fi

if [[ ! -s "$DATA_DIR/pbmc3k_subset.h5ad" ]]; then
    echo "==> Creating a 50% subsample as .h5ad (second pseudo-sample)"
    python3 - "$MTX_DIR" "$DATA_DIR/pbmc3k_subset.h5ad" <<'PY'
import sys
import numpy as np
import scanpy as sc
src, dest = sys.argv[1], sys.argv[2]
a = sc.read_10x_mtx(src, var_names="gene_symbols", cache=False)
a.var_names_make_unique()
rng = np.random.default_rng(0)
idx = np.sort(rng.choice(a.n_obs, size=a.n_obs // 2, replace=False))
a[idx].copy().write_h5ad(dest, compression="gzip")
print(f"wrote {dest}: {len(idx)} cells")
PY
fi

echo "==> Writing sample sheet"
cat > "$DATA_DIR/samples.csv" <<EOF
sample,matrix,condition
PBMC_full,$REPO_ROOT/$MTX_DIR,groupA
PBMC_subset,$REPO_ROOT/$DATA_DIR/pbmc3k_subset.h5ad,groupB
EOF

echo "==> Validating sample sheet before running anything"
python3 scripts/samplesheet.py "$DATA_DIR/samples.csv" matrix

echo "==> Writing test config"
sed -e "s|samplesheet: \"config/samples_example.csv\"|samplesheet: \"$REPO_ROOT/$DATA_DIR/samples.csv\"|" \
    -e "s|outdir: \"results\"|outdir: \"$REPO_ROOT/$DATA_DIR/results\"|" \
    -e 's|project_name: "My single-cell project"|project_name: "scRNA-seq smoke test (PBMC 3k)"|' \
    -e 's|report_pdf: true|report_pdf: false|' \
    -e 's|method: "none"|method: "harmony"|' \
    config/config_local.yaml > "$DATA_DIR/config_test.yaml"

echo "==> Running pipeline"
snakemake -s workflow/Snakefile \
    --configfile "$DATA_DIR/config_test.yaml" \
    --cores "${CORES:-4}" \
    $CONDA_FLAG \
    $DRY_RUN

if [[ -z "$DRY_RUN" ]]; then
    echo
    echo "==> Checking expected outputs"
    fail=0
    for f in "$DATA_DIR/results/h5ad/04_annotated.h5ad" \
             "$DATA_DIR/results/tables/markers.tsv" \
             "$DATA_DIR/results/tables/cluster_annotation.tsv" \
             "$DATA_DIR/results/figures/umap_leiden.png" \
             "$DATA_DIR/results/report/scrnaseq_report.html"; do
        if [[ -s "$f" ]]; then echo "  OK   $f"; else echo "  MISS $f"; fail=1; fi
    done
    for stage in qc doublets cluster markers; do
        f="$DATA_DIR/results/metrics/$stage.json"
        if [[ -s "$f" ]]; then echo "  OK   $f"; else echo "  MISS $f"; fail=1; fi
    done
    echo
    if [[ $fail -eq 0 ]]; then
        n=$(python3 -c "import json;print(json.load(open('$DATA_DIR/results/metrics/cluster.json'))['cluster']['n_clusters'])")
        c=$(python3 -c "import json;print(json.load(open('$DATA_DIR/results/metrics/cluster.json'))['cluster']['cells'])")
        echo "Smoke test PASSED - $c cells, $n clusters, report written."
        echo "Expected cell types for PBMC: T, B, NK, Monocyte, Dendritic, Platelet."
        echo "Labels found:"
        awk -F'\t' 'NR>1{print "  cluster "$1": "$2}' "$DATA_DIR/results/tables/cluster_annotation.tsv"
    else
        echo "Smoke test FAILED - see logs under $DATA_DIR/results/logs/" >&2
        exit 1
    fi
fi
