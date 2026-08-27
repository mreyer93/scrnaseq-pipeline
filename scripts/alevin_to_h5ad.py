"""Convert an alevin-fry / simpleaf quantification directory into an .h5ad.

Kept as a separate script (rather than inline in the rule) so the FASTQ path converges
on exactly the same AnnData interface the matrix path uses, and so it can be run by hand
against an existing alevin-fry output.

Usage: python scripts/alevin_to_h5ad.py <af_quant_dir> <out.h5ad>
"""

import os
import sys


def load_alevin(quant_dir):
    """Read an alevin-fry quant directory into AnnData.

    Prefers pyroe (the tool's own reader, which handles USA-mode splicing layers);
    falls back to reading the matrix-market files directly if pyroe is unavailable.
    """
    try:
        import pyroe
        adata = pyroe.load_fry(quant_dir, output_format="scRNA")
        return adata, "pyroe"
    except ImportError:
        pass
    except Exception as e:
        print(f"pyroe failed to read {quant_dir} ({e}); falling back to direct read")

    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.io import mmread
    from scipy.sparse import csr_matrix

    alevin = os.path.join(quant_dir, "alevin")
    if not os.path.isdir(alevin):
        raise SystemExit(
            f"ERROR: {quant_dir} does not look like an alevin-fry output "
            f"(no 'alevin' subdirectory).")

    mtx = os.path.join(alevin, "quants_mat.mtx")
    rows = os.path.join(alevin, "quants_mat_rows.txt")
    cols = os.path.join(alevin, "quants_mat_cols.txt")
    for p in (mtx, rows, cols):
        if not os.path.exists(p):
            raise SystemExit(f"ERROR: expected file missing from alevin-fry output: {p}")

    X = csr_matrix(mmread(mtx))
    barcodes = pd.read_csv(rows, header=None)[0].astype(str).tolist()
    genes = pd.read_csv(cols, header=None)[0].astype(str).tolist()

    # USA mode emits spliced/unspliced/ambiguous blocks side by side; sum them to get
    # the total per gene, which is what the standard gene-expression workflow expects.
    if len(genes) == 3 * len(set(genes)):
        n = len(set(genes))
        base = genes[:n]
        X = X[:, :n] + X[:, n:2 * n] + X[:, 2 * n:3 * n]
        genes = base

    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes
    adata.var_names = genes
    return adata, "direct"


def main(argv):
    if len(argv) < 3:
        print(__doc__)
        return 2
    quant_dir, out_h5ad = argv[1], argv[2]
    adata, how = load_alevin(quant_dir)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    os.makedirs(os.path.dirname(os.path.abspath(out_h5ad)) or ".", exist_ok=True)
    adata.write_h5ad(out_h5ad, compression="gzip")
    print(f"Read {quant_dir} via {how}: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"Wrote {out_h5ad}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
