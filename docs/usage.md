# Usage

## 1. Install

```bash
mamba env create -f envs/python.yml -n scrnaseq     # analysis (Scanpy)
conda activate scrnaseq
```

Two further environments are used by specific rules and created automatically with
`--use-conda`: `envs/r.yml` (report rendering) and `envs/environment.yml` (simpleaf, only
for `input_mode: fastq`).

## 2. Sample sheet

CSV or TSV; the separator is detected. See the
[README](../README.md#sample-sheets) for accepted column-name aliases.

Matrix mode:

```csv
sample,matrix,condition,batch
CTRL_1,/data/cellranger/CTRL_1/outs/filtered_feature_bc_matrix,control,b1
CTRL_2,/data/CTRL_2.h5,control,b2
TREAT_1,/data/TREAT_1.h5ad,treated,b1
```

FASTQ mode needs both reads — R1 carries the cell barcode and UMI, R2 the cDNA:

```csv
sample,fastq_1,fastq_2,condition
CTRL_1,/data/fq/CTRL_1_R1.fastq.gz,/data/fq/CTRL_1_R2.fastq.gz,control
```

Validate before running anything:

```bash
python3 scripts/samplesheet.py samples.csv matrix     # or: fastq
```

This checks that every input exists and is a recognisable matrix, that no input is shared
between samples, and that sample names are safe for filenames — reporting each problem
with a line number.

## 3. Choosing parameters

The defaults are reasonable for a typical 10x human experiment, but two are worth
deliberate thought:

**QC thresholds.** `qc.min_genes: 200` and `qc.max_pct_mt: 20` are conventional starting
points, not universal truths. Some tissues (heart, muscle) legitimately have high
mitochondrial content, and filtering at 20% would discard real cells. Run once, look at
`figures/qc_before.png`, then set thresholds from what the data actually shows.

**Clustering resolution.** `cluster.resolution` controls how many clusters you get. There
is no correct value. Start at 1.0, then adjust based on whether marker genes distinguish
the clusters you obtained: if two clusters share the same markers, lower it; if one
cluster clearly contains two cell types, raise it.

**Integration.** Leave `integration.method: none` unless the UMAP shows samples separating
where they biologically should not. Harmony removes batch structure, but it can equally
remove real between-sample biology — which matters if sample differences are your
research question.

## 4. Run

```bash
snakemake --configfile my_config.yaml --use-conda --cores 4 -n   # dry run
snakemake --configfile my_config.yaml --use-conda --cores 4
```

Each stage writes its own `.h5ad`, so you can re-run from any point. To retry clustering
at a different resolution without redoing QC and doublets:

```bash
# edit cluster.resolution, then
snakemake --configfile my_config.yaml --use-conda --cores 4 --forcerun cluster
```

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| "All cells were filtered out" | thresholds too strict; check `figures/qc_before.png` |
| "no mitochondrial genes found" | `var_names` are Ensembl IDs, not symbols — the mito filter does nothing; convert IDs or set `max_pct_mt: 100` |
| Only one cluster | raise `cluster.resolution`; or the population really is homogeneous |
| Every cell called a doublet | input is probably raw (unfiltered) droplets rather than called cells |
| Clusters separate by sample, not cell type | batch effect — try `integration.method: harmony` |
| Cell type labels look wrong | marker sets do not match the tissue, or gene symbols do not match the data's naming |
| Out of memory | reduce `cluster.n_hvg`, subsample cells, or use the cloud VM |
| Near-zero barcodes (FASTQ mode) | wrong `chemistry` setting for the kit used |

Each stage writes a log under `<outdir>/logs/` ending in the actual reason for failure
rather than only a traceback.
