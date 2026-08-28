"""Stage 3: normalisation, feature selection, dimensionality reduction, clustering.

Standard Scanpy workflow as described in the single-cell best-practices book:
  shifted-log normalisation -> highly variable genes -> PCA -> (batch integration)
  -> kNN graph -> Leiden clustering -> UMAP

Batch integration with Harmony runs only when there is more than one batch and the user
asked for it. Integration is not free: it can also remove genuine biological differences
between samples, so it is opt-in and the report states whether it was applied.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sc_common import (  # noqa: E402
    StageError, PALETTE, setup_logging, log, save_fig, save_current_fig, save_table,
    save_metrics, ensure_dirs, stage_main,
)


def run():
    import scanpy as sc
    import anndata as ad

    smk = globals()["snakemake"]
    setup_logging(smk.log[0])

    cfg = smk.params["cluster"]
    integration_cfg = smk.params["integration"]
    figures_dir = smk.params["figures_dir"]
    tables_dir = smk.params["tables_dir"]
    metrics_file = smk.params["metrics_file"]
    out_h5ad = smk.output["h5ad"]
    ensure_dirs(figures_dir, tables_dir, os.path.dirname(out_h5ad))

    sc.settings.figdir = figures_dir
    sc.settings.autoshow = False

    adata = ad.read_h5ad(smk.input["h5ad"])
    log(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")
    if adata.n_obs < 30:
        raise StageError(
            f"Only {adata.n_obs} cells remain - too few to cluster meaningfully. "
            f"Check the QC thresholds and the doublet stage.")

    # ------------------------------------------------------- normalisation ---------
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()
    target_sum = cfg.get("target_sum", 1e4)
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.raw = adata          # log-normalised values, used for marker gene plots
    log(f"Normalised to {target_sum:g} counts/cell and log1p-transformed")

    # ---------------------------------------------------- feature selection --------
    n_hvg = int(cfg.get("n_hvg", 2000))
    n_hvg = min(n_hvg, adata.n_vars)
    batch_key = integration_cfg.get("batch_key") or None
    if batch_key and batch_key not in adata.obs.columns:
        raise StageError(
            f"integration.batch_key is {batch_key!r} but that column is not in the data. "
            f"Available: {[c for c in adata.obs.columns if adata.obs[c].dtype.name in ('category','object')]}")
    hvg_batch = batch_key if (batch_key and adata.obs[batch_key].nunique() > 1) else None
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=hvg_batch,
                                    flavor="seurat_v3", layer="counts")
    except Exception as e:
        # seurat_v3 needs raw counts and scikit-misc; fall back to the classic flavor
        log(f"seurat_v3 HVG selection unavailable ({e}); using the default flavor")
        sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, batch_key=hvg_batch)
    n_sel = int(adata.var["highly_variable"].sum())
    log(f"Selected {n_sel} highly variable genes"
        + (f" (per batch: {hvg_batch})" if hvg_batch else ""))
    save_table(
        adata.var.loc[adata.var["highly_variable"]].reset_index()
            .rename(columns={"index": "gene"}).head(500),
        tables_dir, "highly_variable_genes")

    # ------------------------------------------------------------------ PCA --------
    if cfg.get("scale", True):
        sc.pp.scale(adata, max_value=10, zero_center=True)
    n_pcs = int(cfg.get("n_pcs", 50))
    n_pcs = max(2, min(n_pcs, adata.n_obs - 1, n_sel - 1))
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", use_highly_variable=True)
    log(f"PCA: {n_pcs} components")

    fig, ax = plt.subplots(figsize=(6, 3.6))
    vr = adata.uns["pca"]["variance_ratio"]
    ax.plot(range(1, len(vr) + 1), vr * 100, "o-", color=PALETTE[0], ms=3)
    ax.set_xlabel("Principal component")
    ax.set_ylabel("% variance explained")
    ax.set_title("PCA variance", fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_fig(fig, figures_dir, "pca_variance")

    # ---------------------------------------------------------- integration -------
    use_rep = "X_pca"
    integrated = False
    if integration_cfg.get("method", "none") == "harmony" and hvg_batch:
        try:
            # harmonypy is called directly rather than through
            # scanpy.external.pp.harmony_integrate: that wrapper assumes a fixed
            # orientation for harmonypy's corrected matrix, and with current harmonypy
            # releases it writes back a wrongly-shaped array, so integration silently
            # fails. Checking the orientation ourselves works with either convention.
            import harmonypy

            ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, [hvg_batch])
            Z = np.asarray(getattr(ho, "Z_corr"))
            if Z.shape[0] != adata.n_obs:
                Z = Z.T
            if Z.shape[0] != adata.n_obs:
                raise ValueError(
                    f"Harmony returned a {Z.shape} matrix; expected one axis to be "
                    f"{adata.n_obs} cells")
            adata.obsm["X_pca_harmony"] = Z
            use_rep = "X_pca_harmony"
            integrated = True
            log(f"Integrated batches with Harmony on {hvg_batch!r} "
                f"(corrected embedding {Z.shape})")
        except Exception as e:
            log(f"WARNING: Harmony integration failed ({e}); continuing unintegrated. "
                f"Any batch structure in the UMAP is therefore uncorrected.")
    elif integration_cfg.get("method", "none") == "harmony":
        log("Integration requested but there is only one batch - skipping (nothing to integrate)")

    # ------------------------------------------------ neighbours / clustering ------
    n_neighbors = int(cfg.get("n_neighbors", 15))
    n_neighbors = min(n_neighbors, max(2, adata.n_obs - 1))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=use_rep)

    resolution = float(cfg.get("resolution", 1.0))
    try:
        sc.tl.leiden(adata, resolution=resolution, key_added="leiden",
                     flavor="igraph", n_iterations=2, directed=False)
    except TypeError:
        # older scanpy without the igraph flavor argument
        sc.tl.leiden(adata, resolution=resolution, key_added="leiden")
    n_clusters = adata.obs["leiden"].nunique()
    log(f"Leiden clustering at resolution {resolution}: {n_clusters} clusters")
    if n_clusters == 1:
        log("WARNING: only one cluster found. Raise cluster.resolution, or the dataset "
            "may genuinely be homogeneous.")

    sc.tl.umap(adata)
    log("UMAP embedding computed")

    # ------------------------------------------------------------------ plots ------
    _umap_panel(adata, figures_dir, integration_cfg, hvg_batch)

    comp = (adata.obs.groupby(["leiden", "sample"], observed=True).size()
            .reset_index(name="cells"))
    save_table(comp, tables_dir, "cluster_composition")
    sizes = (adata.obs["leiden"].value_counts().rename_axis("cluster")
             .reset_index(name="cells").sort_values("cluster"))
    save_table(sizes, tables_dir, "cluster_sizes")
    _composition_plot(comp, figures_dir)

    adata.write_h5ad(out_h5ad, compression="gzip")
    log(f"Wrote {out_h5ad}")

    save_metrics({"cluster": {
        "n_hvg": n_sel, "n_pcs": n_pcs, "n_neighbors": n_neighbors,
        "resolution": resolution, "n_clusters": int(n_clusters),
        "integration": "harmony" if integrated else "none",
        "batch_key": hvg_batch or "",
        "cells": int(adata.n_obs),
    }}, metrics_file)


def _umap_panel(adata, figures_dir, integration_cfg, batch_key):
    import scanpy as sc
    colour_by = ["leiden", "sample"]
    for extra in ("pct_counts_mt", "n_genes_by_counts", "doublet_score"):
        if extra in adata.obs.columns and adata.obs[extra].notna().any():
            colour_by.append(extra)
    if batch_key and batch_key not in colour_by:
        colour_by.append(batch_key)
    colour_by = [c for c in colour_by if c in adata.obs.columns]

    for c in colour_by:
        try:
            sc.pl.umap(adata, color=c, show=False, frameon=False,
                       legend_loc="on data" if c == "leiden" else "right margin",
                       title=f"UMAP - {c}")
            save_current_fig(figures_dir, f"umap_{c}")
        except Exception as e:
            log(f"  could not plot UMAP coloured by {c}: {e}")
            plt.close("all")


def _composition_plot(comp, figures_dir):
    if comp.empty:
        return
    pivot = comp.pivot(index="leiden", columns="sample", values="cells").fillna(0)
    frac = pivot.div(pivot.sum(axis=1), axis=0)
    fig, ax = plt.subplots(figsize=(9, 4))
    bottom = np.zeros(len(frac))
    for i, s in enumerate(frac.columns):
        ax.bar(frac.index.astype(str), frac[s].values, bottom=bottom,
               label=str(s), color=PALETTE[i % len(PALETTE)])
        bottom += frac[s].values
    ax.set_xlabel("Cluster")
    ax.set_ylabel("Fraction of cells")
    ax.set_title("Sample composition per cluster", fontweight="bold")
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_fig(fig, figures_dir, "cluster_composition")


if __name__ == "__main__":
    stage_main(run)
