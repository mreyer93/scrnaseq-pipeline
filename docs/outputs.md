# Outputs

```
<outdir>/
├── h5ad/
│   ├── 01_qc.h5ad          after cell/gene filtering
│   ├── 02_doublets.h5ad    after doublet detection
│   ├── 03_clustered.h5ad   with PCA, UMAP, Leiden clusters
│   └── 04_annotated.h5ad   with marker results and cell type labels  <- the main object
├── figures/                QC violins, UMAPs, dotplot, heatmaps, composition (PNG)
├── tables/
│   ├── markers.tsv                  top markers per cluster
│   ├── markers_all.tsv              every gene tested per cluster
│   ├── cluster_sizes.tsv            cells per cluster
│   ├── cluster_composition.tsv      cells per cluster per sample
│   ├── cluster_annotation.tsv       provisional cell type labels
│   ├── cell_type_scores.tsv         mean marker-set score per cluster
│   ├── qc_per_sample_before/after.tsv
│   ├── qc_filtering.tsv             cells removed, by reason
│   └── doublets_per_sample.tsv
├── metrics/                one JSON per stage
├── report/                 scrnaseq_report.html / .pdf
└── logs/                   one log per stage
```

## Which file do I actually want?

- **Send to a collaborator/client** → `report/scrnaseq_report.html` (self-contained).
- **Continue the analysis yourself** → `h5ad/04_annotated.h5ad`.
- **Marker genes for a cluster** → `tables/markers.tsv` (top N) or `markers_all.tsv` (all).
- **Cluster proportions** → `tables/cluster_composition.tsv`. See the caveat below before
  comparing these between conditions.

```python
import anndata as ad
adata = ad.read_h5ad("results/h5ad/04_annotated.h5ad")
adata.obs["leiden"]        # cluster assignment
adata.obs["cell_type"]     # provisional label (if marker sets were configured)
adata.obsm["X_umap"]       # UMAP coordinates
adata.layers["counts"]     # raw counts, preserved
```

## markers.tsv columns

| Column | Meaning |
|---|---|
| `cluster` | Leiden cluster the gene marks |
| `gene` | gene symbol |
| `logfoldchange` | log2 fold change, cluster vs all other cells |
| `pval` / `pval_adj` | Wilcoxon p-value and BH-adjusted p-value |
| `score` | test statistic used for ranking |

**These p-values are not valid inference.** Clusters were defined from the same expression
data used to test them ("double dipping"), which inflates significance for every gene.
Use the ranking to identify what distinguishes a cluster; do not report the p-values as
evidence that the cluster itself is real.

## Caveats worth carrying into a write-up

- **Cluster count is a parameter, not a finding.** It follows directly from
  `cluster.resolution`.
- **Cell type labels are provisional** — the best-scoring supplied marker set per cluster,
  flagged *low confidence* when the top two score nearly equally.
- **Cluster proportions between conditions are not tested here.** Apparent compositional
  shifts are frequently noise with few biological replicates; use a compositional method
  (scCODA, miloR) if that comparison matters.
- **Ambient RNA is not corrected**, so a low level of cross-cluster "expression" of highly
  expressed genes is expected.
