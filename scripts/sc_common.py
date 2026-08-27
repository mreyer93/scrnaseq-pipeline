"""Shared helpers for the single-cell analysis stages.

Each stage script reads an .h5ad, does one conceptual step, and writes an .h5ad plus any
figures (PNG), tables (TSV) and metrics (JSON) it produced. Keeping the stages separate
means an expensive run can be resumed part-way, and each step's output can be inspected
on its own rather than only at the end.
"""

import json
import os
import sys
import traceback

import matplotlib
matplotlib.use("Agg")            # no display on servers/CI
import matplotlib.pyplot as plt  # noqa: E402


PALETTE = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
           "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD"]


class StageError(RuntimeError):
    """Raised with a message intended to be read directly by the user."""


def setup_logging(log_path):
    """Tee stdout/stderr to a log file so failures are diagnosable after the fact."""
    if not log_path:
        return None
    os.makedirs(os.path.dirname(os.path.abspath(log_path)) or ".", exist_ok=True)
    fh = open(log_path, "w", encoding="utf-8")

    class _Tee:
        def __init__(self, *streams):
            self.streams = streams

        def write(self, data):
            for s in self.streams:
                s.write(data)
                s.flush()

        def flush(self):
            for s in self.streams:
                s.flush()

    sys.stdout = _Tee(sys.__stdout__, fh)
    sys.stderr = _Tee(sys.__stderr__, fh)
    return fh


def log(msg):
    print(msg, flush=True)


def die(msg):
    """Fail with a message the user can act on, not a bare traceback."""
    print(f"\nERROR: {msg}\n", file=sys.stderr)
    raise SystemExit(1)


def ensure_dirs(*paths):
    for p in paths:
        os.makedirs(p, exist_ok=True)


def save_fig(fig, figures_dir, name, dpi=150):
    ensure_dirs(figures_dir)
    path = os.path.join(figures_dir, f"{name}.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    log(f"  figure: {path}")
    return path


def save_current_fig(figures_dir, name, dpi=150):
    """Save whatever scanpy just drew (scanpy plots onto the current figure)."""
    fig = plt.gcf()
    return save_fig(fig, figures_dir, name, dpi=dpi)


def save_table(df, tables_dir, name, index=False):
    ensure_dirs(tables_dir)
    path = os.path.join(tables_dir, f"{name}.tsv")
    df.to_csv(path, sep="\t", index=index)
    log(f"  table: {path} ({len(df)} rows)")
    return path


def save_metrics(metrics, path):
    """Write this stage's metrics to its own file.

    One file per stage, rather than all stages appending to a shared metrics.json:
    Snakemake deletes a rule's declared outputs before rerunning it, so a shared file
    would silently lose the earlier stages' sections whenever a later stage reran.
    The report merges the per-stage files at the end.
    """
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(metrics, fh, indent=2, default=str)
    log(f"  metrics: {path}")


def load_matrix(path, kind, sample_name):
    """Read one sample's count matrix into AnnData, whatever form it arrived in."""
    import scanpy as sc
    import anndata as ad

    log(f"Loading {sample_name} ({kind}): {path}")
    try:
        if kind == "10x_mtx":
            adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
        elif kind == "h5":
            adata = sc.read_10x_h5(path)
        elif kind == "h5ad":
            adata = ad.read_h5ad(path)
        else:
            raise StageError(f"unknown matrix kind {kind!r} for {path}")
    except Exception as e:
        raise StageError(
            f"Failed to read {sample_name} from {path} (detected as {kind}): {e}\n"
            f"For a 10x directory the folder must contain matrix.mtx(.gz), "
            f"barcodes.tsv(.gz) and features.tsv(.gz) or genes.tsv(.gz)."
        )

    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise StageError(f"{sample_name}: matrix is empty ({adata.shape})")
    log(f"  {adata.n_obs} cells x {adata.n_vars} genes")
    return adata


def guess_species_prefixes(var_names):
    """Mitochondrial/ribosomal gene prefixes differ by species and annotation.

    Human symbols are upper case (MT-, RPS/RPL), mouse are title case (mt-, Rps/Rpl),
    and Ensembl IDs carry no prefix at all. Detect rather than assume, so QC does not
    silently compute 0% mitochondrial content on a mouse dataset.
    """
    names = list(var_names[: min(len(var_names), 5000)])
    upper = sum(1 for n in names if n.isupper())
    has_mt_upper = any(n.startswith("MT-") for n in names)
    has_mt_lower = any(n.startswith("mt-") for n in names)
    if has_mt_upper:
        return {"mito": ("MT-",), "ribo": ("RPS", "RPL"), "hb": ("HBA", "HBB")}
    if has_mt_lower:
        return {"mito": ("mt-",), "ribo": ("Rps", "Rpl"), "hb": ("Hba", "Hbb")}
    # Fall back on case heuristics; may still find nothing (e.g. Ensembl IDs)
    if upper > len(names) / 2:
        return {"mito": ("MT-", "MT."), "ribo": ("RPS", "RPL"), "hb": ("HBA", "HBB")}
    return {"mito": ("mt-", "Mt-"), "ribo": ("Rps", "Rpl"), "hb": ("Hba", "Hbb")}


def stage_main(fn):
    """Run a stage, converting expected failures into clean messages.

    Snakemake surfaces the rule name and log path on failure; what we add here is that
    the log's last lines say what actually went wrong rather than ending in a traceback
    through library internals.
    """
    try:
        fn()
    except StageError as e:
        die(str(e))
    except MemoryError:
        die("Ran out of memory.\n"
            "Single-cell objects are large. Options: raise the machine's RAM, reduce the "
            "dataset (subsample cells), or run the cloud config on a bigger VM.")
    except Exception:
        traceback.print_exc()
        die("Unexpected failure - the traceback above is the primary evidence. "
            "If it mentions a missing gene/column, check that the sample sheet metadata "
            "matches the data.")
