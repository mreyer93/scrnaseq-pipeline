"""Stage 1: load every sample, compute QC metrics, filter cells and genes.

Follows the quality-control chapter of the single-cell best-practices book
(https://www.sc-best-practices.org/): per-cell counts, detected genes, and
mitochondrial fraction, filtered with thresholds that are explicit and configurable
rather than silently hard-coded.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sc_common import (  # noqa: E402
    StageError, PALETTE, setup_logging, log, save_fig, save_table, save_metrics,
    ensure_dirs, load_matrix, guess_species_prefixes, stage_main,
)


def run():
    import scanpy as sc
    import anndata as ad

    smk = globals()["snakemake"]
    setup_logging(smk.log[0])

    samples = smk.params["samples"]          # dict: name -> {matrix, kind, metadata...}
    extra_cols = smk.params["extra_cols"]
    qc_cfg = smk.params["qc"]
    figures_dir = smk.params["figures_dir"]
    tables_dir = smk.params["tables_dir"]
    metrics_file = smk.params["metrics_file"]
    out_h5ad = smk.output["h5ad"]

    ensure_dirs(figures_dir, tables_dir, os.path.dirname(out_h5ad))

    # ---------------------------------------------------------------- load ---------
    adatas = {}
    for name, s in samples.items():
        a = load_matrix(s["matrix"], s.get("kind", "10x_mtx"), name)
        for c in extra_cols:
            a.obs[c] = str(s.get(c, ""))
        adatas[name] = a

    if len(adatas) == 1:
        adata = next(iter(adatas.values()))
        adata.obs["sample"] = next(iter(adatas))
    else:
        # Concatenating on the intersection of genes: samples processed with different
        # reference versions otherwise produce a matrix full of NaNs.
        adata = ad.concat(adatas, label="sample", join="inner", index_unique="-")
        n_genes = {k: v.n_vars for k, v in adatas.items()}
        log(f"Concatenated {len(adatas)} samples on shared genes: "
            f"{adata.n_vars} kept (per-sample: {n_genes})")
        if adata.n_vars < 0.5 * max(n_genes.values()):
            log("WARNING: fewer than half the genes are shared across samples. "
                "This usually means they were quantified against different references.")
    adata.obs["sample"] = adata.obs["sample"].astype("category")
    log(f"Combined: {adata.n_obs} cells x {adata.n_vars} genes")

    if adata.n_obs == 0:
        raise StageError("No cells loaded.")

    # ------------------------------------------------------- QC metrics ------------
    prefixes = guess_species_prefixes(adata.var_names)
    adata.var["mt"] = adata.var_names.str.startswith(prefixes["mito"])
    adata.var["ribo"] = adata.var_names.str.startswith(prefixes["ribo"])
    adata.var["hb"] = adata.var_names.str.startswith(prefixes["hb"])
    n_mt = int(adata.var["mt"].sum())
    log(f"Mitochondrial genes detected: {n_mt} (prefixes {prefixes['mito']})")
    if n_mt == 0:
        log("WARNING: no mitochondrial genes found. The mitochondrial-fraction filter "
            "will not remove anything. This is expected if var_names are Ensembl IDs "
            "rather than gene symbols.")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"],
                               inplace=True, percent_top=None, log1p=False)

    pre = pd.DataFrame({
        "sample": adata.obs["sample"].values,
        "n_genes": adata.obs["n_genes_by_counts"].values,
        "total_counts": adata.obs["total_counts"].values,
        "pct_mt": adata.obs["pct_counts_mt"].values,
    })
    per_sample = pre.groupby("sample", observed=True).agg(
        cells=("n_genes", "size"),
        median_genes=("n_genes", "median"),
        median_counts=("total_counts", "median"),
        median_pct_mt=("pct_mt", "median"),
    ).reset_index()
    save_table(per_sample, tables_dir, "qc_per_sample_before")

    # ------------------------------------------------------------- figures ---------
    _qc_violin(adata, figures_dir, "qc_before")

    # ------------------------------------------------------------ filtering --------
    min_genes = int(qc_cfg.get("min_genes", 200))
    max_genes = qc_cfg.get("max_genes", None)
    min_counts = int(qc_cfg.get("min_counts", 0))
    max_pct_mt = float(qc_cfg.get("max_pct_mt", 20))
    min_cells = int(qc_cfg.get("min_cells_per_gene", 3))

    n0 = adata.n_obs
    reasons = {}

    keep = adata.obs["n_genes_by_counts"] >= min_genes
    reasons[f"n_genes < {min_genes}"] = int((~keep).sum())
    adata = adata[keep].copy()

    if min_counts > 0:
        keep = adata.obs["total_counts"] >= min_counts
        reasons[f"total_counts < {min_counts}"] = int((~keep).sum())
        adata = adata[keep].copy()

    if max_genes:
        keep = adata.obs["n_genes_by_counts"] <= int(max_genes)
        reasons[f"n_genes > {max_genes}"] = int((~keep).sum())
        adata = adata[keep].copy()

    if n_mt > 0:
        keep = adata.obs["pct_counts_mt"] <= max_pct_mt
        reasons[f"pct_mt > {max_pct_mt}"] = int((~keep).sum())
        adata = adata[keep].copy()

    n_genes_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    genes_removed = n_genes_before - adata.n_vars

    log(f"Filtering removed {n0 - adata.n_obs} / {n0} cells:")
    for k, v in reasons.items():
        log(f"  {v:>7} cells: {k}")
    log(f"  {genes_removed} genes detected in fewer than {min_cells} cells")
    log(f"Remaining: {adata.n_obs} cells x {adata.n_vars} genes")

    if adata.n_obs == 0:
        raise StageError(
            f"All {n0} cells were filtered out. The thresholds are too strict for this "
            f"dataset (removals: {reasons}). Inspect figures/qc_before.png and relax "
            f"qc.min_genes / qc.max_pct_mt in the config.")
    if adata.n_obs < 50:
        log(f"WARNING: only {adata.n_obs} cells survived filtering. Clustering results "
            f"from so few cells are unlikely to be meaningful.")

    _qc_violin(adata, figures_dir, "qc_after")

    post = pd.DataFrame({
        "sample": adata.obs["sample"].values,
        "n_genes": adata.obs["n_genes_by_counts"].values,
        "total_counts": adata.obs["total_counts"].values,
        "pct_mt": adata.obs["pct_counts_mt"].values,
    })
    per_sample_after = post.groupby("sample", observed=True).agg(
        cells=("n_genes", "size"),
        median_genes=("n_genes", "median"),
        median_counts=("total_counts", "median"),
        median_pct_mt=("pct_mt", "median"),
    ).reset_index()
    save_table(per_sample_after, tables_dir, "qc_per_sample_after")

    filt = pd.DataFrame({"reason": list(reasons) + [f"gene in <{min_cells} cells"],
                         "removed": list(reasons.values()) + [genes_removed]})
    save_table(filt, tables_dir, "qc_filtering")

    # keep the raw counts around; later stages normalise in place
    adata.layers["counts"] = adata.X.copy()
    adata.write_h5ad(out_h5ad, compression="gzip")
    log(f"Wrote {out_h5ad}")

    save_metrics({
        "qc": {
            "cells_before": int(n0),
            "cells_after": int(adata.n_obs),
            "genes_after": int(adata.n_vars),
            "n_samples": int(adata.obs['sample'].nunique()),
            "mito_genes_found": n_mt,
            "thresholds": {"min_genes": min_genes, "max_genes": max_genes,
                           "min_counts": min_counts, "max_pct_mt": max_pct_mt,
                           "min_cells_per_gene": min_cells},
        }
    }, metrics_file)


def _qc_violin(adata, figures_dir, name):
    """Violin plots of the three standard QC metrics, split by sample."""
    metrics = [("n_genes_by_counts", "Genes per cell"),
               ("total_counts", "UMI counts per cell"),
               ("pct_counts_mt", "% mitochondrial")]
    samples = list(adata.obs["sample"].cat.categories) if hasattr(
        adata.obs["sample"], "cat") else sorted(adata.obs["sample"].unique())
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, (col, title) in zip(axes, metrics):
        data = [adata.obs.loc[adata.obs["sample"] == s, col].values for s in samples]
        data = [d for d in data if len(d) > 0]
        if not data:
            ax.set_visible(False)
            continue
        parts = ax.violinplot(data, showmedians=True, widths=0.85)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(PALETTE[i % len(PALETTE)])
            body.set_alpha(0.75)
        ax.set_xticks(range(1, len(data) + 1))
        ax.set_xticklabels(samples[:len(data)], rotation=45, ha="right", fontsize=8)
        ax.set_title(title, fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_fig(fig, figures_dir, name)


if __name__ == "__main__":
    stage_main(run)
