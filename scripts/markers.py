"""Stage 4: marker genes per cluster, and optional cell-type scoring.

Marker detection uses the Wilcoxon rank-sum test, Scanpy's default and the usual choice
for droplet data. A caveat worth stating plainly in any write-up: these p-values are
computed on the same cells that were used to define the clusters, so they are inflated
and should be read as a ranking, not as evidence that a cluster is real.

If `annotation.marker_sets` is supplied, each cluster is also scored against those gene
sets (scanpy's score_genes) and given a best-match label. That is a starting hypothesis
for annotation, not an answer - it is only ever as good as the supplied marker sets.
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

    cfg = smk.params["markers"]
    annotation_cfg = smk.params["annotation"] or {}
    figures_dir = smk.params["figures_dir"]
    tables_dir = smk.params["tables_dir"]
    metrics_file = smk.params["metrics_file"]
    out_h5ad = smk.output["h5ad"]
    ensure_dirs(figures_dir, tables_dir, os.path.dirname(out_h5ad))

    sc.settings.figdir = figures_dir
    sc.settings.autoshow = False

    adata = ad.read_h5ad(smk.input["h5ad"])
    log(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    if "leiden" not in adata.obs.columns:
        raise StageError("No 'leiden' clustering found - the clustering stage must run first.")

    n_clusters = adata.obs["leiden"].nunique()
    if n_clusters < 2:
        log("Only one cluster - marker detection needs at least two groups to compare. "
            "Writing an empty marker table and continuing.")
        save_table(pd.DataFrame(columns=["cluster", "gene", "logfoldchange",
                                         "pval", "pval_adj", "score"]),
                   tables_dir, "markers")
        adata.write_h5ad(out_h5ad, compression="gzip")
        save_metrics({"markers": {"n_clusters": int(n_clusters), "n_markers": 0}},
                     metrics_file)
        return

    method = cfg.get("method", "wilcoxon")
    n_top = int(cfg.get("n_top", 25))

    # Marker tests run on log-normalised values, not the scaled matrix used for PCA
    sc.tl.rank_genes_groups(adata, groupby="leiden", method=method, use_raw=True)
    log(f"Ranked marker genes per cluster ({method})")

    df = sc.get.rank_genes_groups_df(adata, group=None)
    df = df.rename(columns={"group": "cluster", "names": "gene",
                            "logfoldchanges": "logfoldchange", "scores": "score"})
    keep = ["cluster", "gene", "logfoldchange", "pvals", "pvals_adj", "score"]
    df = df[[c for c in keep if c in df.columns]].rename(
        columns={"pvals": "pval", "pvals_adj": "pval_adj"})
    save_table(df, tables_dir, "markers_all")

    padj_cut = float(cfg.get("padj", 0.05))
    lfc_cut = float(cfg.get("min_logfc", 0.25))
    sig = df[(df["pval_adj"] < padj_cut) & (df["logfoldchange"] >= lfc_cut)]
    top = (sig.sort_values(["cluster", "pval_adj"])
              .groupby("cluster", observed=True).head(n_top))
    save_table(top, tables_dir, "markers")
    log(f"{len(sig)} significant markers (padj<{padj_cut}, logFC>={lfc_cut}); "
        f"top {n_top} per cluster saved")

    # ------------------------------------------------------------------ plots ------
    top_genes = (top.groupby("cluster", observed=True).head(5)["gene"]
                 .drop_duplicates().tolist())
    top_genes = [g for g in top_genes if g in (adata.raw.var_names if adata.raw
                                               else adata.var_names)]
    if top_genes:
        try:
            sc.pl.dotplot(adata, var_names=top_genes[:40], groupby="leiden",
                          show=False, standard_scale="var")
            save_current_fig(figures_dir, "marker_dotplot")
        except Exception as e:
            log(f"  dotplot failed: {e}")
            plt.close("all")
        try:
            sc.pl.rank_genes_groups_heatmap(adata, n_genes=4, show=False,
                                            show_gene_labels=True, standard_scale="var")
            save_current_fig(figures_dir, "marker_heatmap")
        except Exception as e:
            log(f"  marker heatmap failed: {e}")
            plt.close("all")
        for g in top_genes[:6]:
            try:
                sc.pl.umap(adata, color=g, show=False, frameon=False, cmap="viridis")
                save_current_fig(figures_dir, f"umap_gene_{_safe(g)}")
            except Exception:
                plt.close("all")

    # ------------------------------------------------------- cell type scoring -----
    marker_sets = annotation_cfg.get("marker_sets") or {}
    labels = None
    if marker_sets:
        labels = _score_cell_types(adata, marker_sets, tables_dir, figures_dir)
        if labels is not None:
            adata.obs["cell_type"] = adata.obs["leiden"].map(labels).astype("category")
            try:
                import scanpy as sc2
                sc2.pl.umap(adata, color="cell_type", show=False, frameon=False,
                            title="UMAP - predicted cell type")
                save_current_fig(figures_dir, "umap_cell_type")
            except Exception:
                plt.close("all")
            save_table(pd.DataFrame({"cluster": list(labels), "cell_type": list(labels.values())}),
                       tables_dir, "cluster_annotation")
    else:
        log("No annotation.marker_sets configured - clusters are left unlabelled. "
            "Supply marker sets in the config to get provisional labels.")

    adata.write_h5ad(out_h5ad, compression="gzip")
    log(f"Wrote {out_h5ad}")

    save_metrics({"markers": {
        "method": method,
        "n_clusters": int(n_clusters),
        "n_markers": int(len(sig)),
        "padj": padj_cut,
        "min_logfc": lfc_cut,
        "annotated": bool(labels),
    }}, metrics_file)


def _score_cell_types(adata, marker_sets, tables_dir, figures_dir):
    """Score each cluster against user-supplied marker gene sets."""
    import scanpy as sc

    available = set(adata.raw.var_names if adata.raw is not None else adata.var_names)
    usable = {}
    for label, genes in marker_sets.items():
        present = [g for g in genes if g in available]
        if not present:
            log(f"  marker set {label!r}: none of its genes are in the data - skipped")
            continue
        if len(present) < len(genes):
            log(f"  marker set {label!r}: using {len(present)}/{len(genes)} genes "
                f"(missing: {sorted(set(genes) - set(present))[:5]})")
        usable[label] = present

    if not usable:
        log("None of the configured marker sets matched any genes - skipping annotation. "
            "Check that the marker gene symbols match the dataset's gene naming.")
        return None

    for label, genes in usable.items():
        sc.tl.score_genes(adata, gene_list=genes, score_name=f"score_{label}")

    score_cols = [f"score_{l}" for l in usable]
    per_cluster = adata.obs.groupby("leiden", observed=True)[score_cols].mean()
    per_cluster.columns = [c.replace("score_", "") for c in per_cluster.columns]
    save_table(per_cluster.reset_index(), tables_dir, "cell_type_scores")

    best = per_cluster.idxmax(axis=1)
    margin = per_cluster.max(axis=1) - per_cluster.apply(
        lambda r: r.nlargest(2).iloc[-1] if len(r) > 1 else r.max(), axis=1)
    labels = {}
    for cl in per_cluster.index:
        lab = best[cl]
        # a near-tie between two sets is not a confident call; say so rather than
        # presenting an arbitrary winner as an answer
        labels[cl] = lab if margin[cl] > 0.01 else f"{lab} (low confidence)"
    log("Provisional cluster labels: " + ", ".join(f"{k}={v}" for k, v in labels.items()))

    fig, ax = plt.subplots(figsize=(max(6, 0.8 * len(per_cluster.columns) + 3),
                                    max(3, 0.35 * len(per_cluster) + 2)))
    im = ax.imshow(per_cluster.values, aspect="auto", cmap="RdBu_r")
    ax.set_xticks(range(len(per_cluster.columns)))
    ax.set_xticklabels(per_cluster.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(per_cluster.index)))
    ax.set_yticklabels(per_cluster.index, fontsize=8)
    ax.set_xlabel("Marker set")
    ax.set_ylabel("Cluster")
    ax.set_title("Mean marker-set score per cluster", fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.8, label="score")
    fig.tight_layout()
    save_fig(fig, figures_dir, "cell_type_scores")
    return labels


def _safe(name):
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in str(name))


if __name__ == "__main__":
    stage_main(run)
