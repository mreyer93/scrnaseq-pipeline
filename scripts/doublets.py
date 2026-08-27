"""Stage 2: doublet detection with Scrublet.

Droplet protocols co-encapsulate two cells often enough that undetected doublets show up
as spurious "intermediate" clusters that look like novel transitional cell states. Scrublet
simulates doublets from the observed data and scores each barcode against them.

Scoring is done per sample: doublet rate scales with loading concentration, so pooling
samples before scoring would apply one threshold to libraries with genuinely different
rates. Flagged cells are removed by default, but `doublets.remove: false` keeps them and
only annotates - useful when you want to see where they landed before committing.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sc_common import (  # noqa: E402
    StageError, PALETTE, setup_logging, log, save_fig, save_table, save_metrics,
    ensure_dirs, stage_main,
)


def run():
    import scanpy as sc
    import anndata as ad

    smk = globals()["snakemake"]
    setup_logging(smk.log[0])

    cfg = smk.params["doublets"]
    figures_dir = smk.params["figures_dir"]
    tables_dir = smk.params["tables_dir"]
    metrics_file = smk.params["metrics_file"]
    out_h5ad = smk.output["h5ad"]
    ensure_dirs(figures_dir, tables_dir, os.path.dirname(out_h5ad))

    adata = ad.read_h5ad(smk.input["h5ad"])
    log(f"Loaded {adata.n_obs} cells x {adata.n_vars} genes")

    if not cfg.get("enabled", True):
        log("Doublet detection disabled in config - passing data through unchanged.")
        adata.obs["predicted_doublet"] = False
        adata.obs["doublet_score"] = np.nan
        adata.write_h5ad(out_h5ad, compression="gzip")
        save_metrics({"doublets": {"enabled": False}}, metrics_file)
        return

    expected_rate = float(cfg.get("expected_rate", 0.06))
    threshold = cfg.get("threshold", None)
    remove = bool(cfg.get("remove", True))

    scores = np.full(adata.n_obs, np.nan)
    predicted = np.zeros(adata.n_obs, dtype=bool)
    per_sample_rows = []

    for sample in adata.obs["sample"].unique():
        mask = (adata.obs["sample"] == sample).values
        sub = adata[mask].copy()
        n = sub.n_obs
        if n < 30:
            log(f"{sample}: only {n} cells - too few for Scrublet, skipping "
                f"(simulation needs a reasonable population to model)")
            per_sample_rows.append({"sample": sample, "cells": n, "doublets": 0,
                                    "pct_doublets": 0.0, "note": "skipped (too few cells)"})
            continue
        log(f"{sample}: scoring {n} cells (expected rate {expected_rate})")
        try:
            # scanpy's wrapper writes predicted_doublet / doublet_score into sub.obs
            sc.pp.scrublet(sub, expected_doublet_rate=expected_rate,
                           threshold=threshold, verbose=False)
        except Exception as e:
            log(f"  Scrublet failed for {sample}: {e}")
            log("  continuing without doublet calls for this sample")
            per_sample_rows.append({"sample": sample, "cells": n, "doublets": 0,
                                    "pct_doublets": 0.0, "note": f"failed: {e}"})
            continue

        s = sub.obs["doublet_score"].to_numpy()
        p = sub.obs["predicted_doublet"].to_numpy().astype(bool)
        scores[mask] = s
        predicted[mask] = p
        pct = 100.0 * p.sum() / n
        log(f"  {int(p.sum())} doublets ({pct:.1f}%)")
        per_sample_rows.append({"sample": sample, "cells": n, "doublets": int(p.sum()),
                                "pct_doublets": round(pct, 2), "note": ""})

    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = predicted

    tbl = pd.DataFrame(per_sample_rows)
    save_table(tbl, tables_dir, "doublets_per_sample")

    _plot_scores(adata, figures_dir)

    n_before = adata.n_obs
    if remove:
        adata = adata[~adata.obs["predicted_doublet"]].copy()
        log(f"Removed {n_before - adata.n_obs} predicted doublets; "
            f"{adata.n_obs} cells remain")
        if adata.n_obs == 0:
            raise StageError(
                "Every cell was called a doublet. This is not a real result - it usually "
                "means the input is not droplet-based single-cell data, or the matrix is "
                "unfiltered raw droplets. Set doublets.enabled: false to bypass.")
    else:
        log(f"Doublets annotated but not removed ({int(predicted.sum())} flagged)")

    adata.write_h5ad(out_h5ad, compression="gzip")
    log(f"Wrote {out_h5ad}")

    save_metrics({"doublets": {
        "enabled": True,
        "expected_rate": expected_rate,
        "flagged": int(predicted.sum()),
        "removed": int(n_before - adata.n_obs) if remove else 0,
        "cells_after": int(adata.n_obs),
    }}, metrics_file)


def _plot_scores(adata, figures_dir):
    samples = sorted(adata.obs["sample"].unique())
    fig, ax = plt.subplots(figsize=(9, 4))
    plotted = False
    for i, s in enumerate(samples):
        v = adata.obs.loc[adata.obs["sample"] == s, "doublet_score"].dropna().values
        if len(v) == 0:
            continue
        ax.hist(v, bins=50, alpha=0.55, label=str(s), color=PALETTE[i % len(PALETTE)])
        plotted = True
    if not plotted:
        plt.close(fig)
        return
    ax.set_xlabel("Scrublet doublet score")
    ax.set_ylabel("Cells")
    ax.set_title("Doublet score distribution", fontweight="bold")
    ax.legend(fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save_fig(fig, figures_dir, "doublet_scores")


if __name__ == "__main__":
    stage_main(run)
