"""Shared setup: config defaults, sample sheet loading, pre-flight validation."""

import os
import sys
from os.path import join

WORKFLOW_DIR = str(workflow.basedir)
REPO_ROOT = os.path.dirname(WORKFLOW_DIR)
SCRIPTS_DIR = join(REPO_ROOT, "scripts")
sys.path.insert(0, SCRIPTS_DIR)

from samplesheet import (  # noqa: E402
    load_samplesheet, summarise, SampleSheetError,
)

# ---------------------------------------------------------------- config defaults ---
config.setdefault("outdir", "results")
config.setdefault("input_mode", "matrix")        # matrix | fastq
config.setdefault("project_name", "Single-cell RNA-seq analysis")
config.setdefault("make_report", True)
config.setdefault("report_pdf", True)

config.setdefault("qc", {})
config["qc"].setdefault("min_genes", 200)
config["qc"].setdefault("max_genes", None)
config["qc"].setdefault("min_counts", 0)
config["qc"].setdefault("max_pct_mt", 20)
config["qc"].setdefault("min_cells_per_gene", 3)

config.setdefault("doublets", {})
config["doublets"].setdefault("enabled", True)
config["doublets"].setdefault("expected_rate", 0.06)
config["doublets"].setdefault("threshold", None)
config["doublets"].setdefault("remove", True)

config.setdefault("cluster", {})
config["cluster"].setdefault("target_sum", 1e4)
config["cluster"].setdefault("n_hvg", 2000)
config["cluster"].setdefault("scale", True)
config["cluster"].setdefault("n_pcs", 50)
config["cluster"].setdefault("n_neighbors", 15)
config["cluster"].setdefault("resolution", 1.0)

config.setdefault("integration", {})
config["integration"].setdefault("method", "none")   # none | harmony
config["integration"].setdefault("batch_key", "sample")

config.setdefault("markers", {})
config["markers"].setdefault("method", "wilcoxon")
config["markers"].setdefault("n_top", 25)
config["markers"].setdefault("padj", 0.05)
config["markers"].setdefault("min_logfc", 0.25)

config.setdefault("annotation", {})
config.setdefault("threads", {})
config["threads"].setdefault("quant", 8)

OUTDIR = config["outdir"]
INPUT_MODE = config["input_mode"]
if INPUT_MODE not in ("matrix", "fastq"):
    raise WorkflowError(f"config 'input_mode' must be 'matrix' or 'fastq', got {INPUT_MODE!r}")

FIGURES_DIR = join(OUTDIR, "figures")
TABLES_DIR = join(OUTDIR, "tables")
METRICS_DIR = join(OUTDIR, "metrics")
def metrics_path(stage):
    return join(METRICS_DIR, f"{stage}.json")

# ------------------------------------------------------------------- sample sheet ---
if "samplesheet" not in config:
    raise WorkflowError(
        "config is missing 'samplesheet'. It needs at least a sample-name column and "
        "either a 'matrix' column (input_mode: matrix) or fastq_1/fastq_2 "
        "(input_mode: fastq). See config/samples_example.csv.")

try:
    SAMPLES, EXTRA_COLS = load_samplesheet(config["samplesheet"],
                                           input_mode=INPUT_MODE, check_files=True)
except SampleSheetError as e:
    raise WorkflowError(str(e))

SAMPLE_NAMES = list(SAMPLES)

# Integration only makes sense across more than one batch
_batch_key = config["integration"].get("batch_key") or "sample"
if _batch_key == "sample":
    _n_batches = len(SAMPLE_NAMES)
elif _batch_key in EXTRA_COLS:
    _n_batches = len({SAMPLES[s].get(_batch_key, "") for s in SAMPLE_NAMES})
else:
    raise WorkflowError(
        f"integration.batch_key is {_batch_key!r}, which is neither 'sample' nor a "
        f"column in the sample sheet (columns: {EXTRA_COLS})")

if config["integration"]["method"] == "harmony" and _n_batches < 2:
    print(f"NOTE: integration.method is 'harmony' but {_batch_key!r} has only "
          f"{_n_batches} level - integration will be skipped (nothing to integrate).")


onstart:
    print("=" * 72)
    print(f"  {config['project_name']}")
    print("=" * 72)
    print(summarise(SAMPLES, EXTRA_COLS, INPUT_MODE))
    print(f"  output dir  : {OUTDIR}")
    print(f"  QC          : min_genes={config['qc']['min_genes']}, "
          f"max_pct_mt={config['qc']['max_pct_mt']}")
    print(f"  doublets    : {'on' if config['doublets']['enabled'] else 'off'}"
          f"{' (annotate only)' if not config['doublets']['remove'] else ''}")
    print(f"  clustering  : resolution={config['cluster']['resolution']}, "
          f"n_hvg={config['cluster']['n_hvg']}")
    print(f"  integration : {config['integration']['method']} on {_batch_key!r} "
          f"({_n_batches} batches)")
    print("=" * 72)


def final_outputs():
    out = [
        join(OUTDIR, "h5ad", "04_annotated.h5ad"),
        join(TABLES_DIR, "markers.tsv"),
        metrics_path("markers"),
    ]
    if config["make_report"]:
        out.append(join(OUTDIR, "report", "scrnaseq_report.html"))
        if config["report_pdf"]:
            out.append(join(OUTDIR, "report", "scrnaseq_report.pdf"))
    return out
