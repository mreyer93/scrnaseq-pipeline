"""Sample sheet parsing and validation for the single-cell RNA-seq pipeline.

Two input modes, chosen per run by config['input_mode']:

  matrix   Each sample points at an existing count matrix. This is what you usually get
           from a sequencing core or a public dataset, and it is the local/laptop path.
           Accepted forms, auto-detected:
             * a 10x Cell Ranger directory containing matrix.mtx(.gz) + barcodes + features
             * an .h5 file (Cell Ranger filtered_feature_bc_matrix.h5)
             * an .h5ad file (AnnData)

  fastq    Each sample points at raw FASTQs, quantified with simpleaf/alevin-fry.
           This is the cloud/server path.

Usable two ways:
  * imported by the workflow (`load_samplesheet`)
  * run directly as a pre-flight check:
        python scripts/samplesheet.py samples.csv matrix
"""

import csv
import os
import re
import sys
from collections import OrderedDict, defaultdict

_ALIASES = {
    "sample":    ["sample", "sampleid", "samplename", "name", "id", "library"],
    "matrix":    ["matrix", "path", "matrixpath", "counts", "countmatrix", "data",
                  "directory", "dir", "h5", "h5ad", "file"],
    "fastq_1":   ["fastq1", "fastqr1", "r1", "read1", "reads1", "fq1", "file1"],
    "fastq_2":   ["fastq2", "fastqr2", "r2", "read2", "reads2", "fq2", "file2"],
    "condition": ["condition", "group", "treatment", "class", "phenotype", "status",
                  "genotype", "timepoint"],
    "batch":     ["batch", "run", "lane", "flowcell", "plate", "donor", "patient",
                  "subject"],
}

_SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
_TENX_MTX = ("matrix.mtx", "matrix.mtx.gz")


class SampleSheetError(ValueError):
    """Raised with a message intended to be read directly by the user."""


def _norm_key(k):
    return re.sub(r"[^a-z0-9]", "", (k or "").strip().lower())


def _sniff_delimiter(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        head = fh.readline()
    if not head:
        raise SampleSheetError(f"{path} is empty.")
    counts = {d: head.count(d) for d in ("\t", ",", ";")}
    delim = max(counts, key=counts.get)
    if counts[delim] == 0:
        raise SampleSheetError(
            f"Could not find a tab, comma or semicolon separator in the header of {path}.\n"
            f"Header was: {head.strip()!r}")
    return delim


def _map_columns(fieldnames):
    mapping, unmapped = {}, []
    for raw in fieldnames:
        key = _norm_key(raw)
        for canonical, aliases in _ALIASES.items():
            if key == _norm_key(canonical) or key in aliases:
                if canonical not in mapping.values():
                    mapping[raw] = canonical
                    break
        else:
            unmapped.append(raw)
    return mapping, unmapped


def detect_matrix_kind(path):
    """Return '10x_mtx' | 'h5' | 'h5ad' | None for a given matrix path."""
    if os.path.isdir(path):
        names = set(os.listdir(path))
        if any(m in names for m in _TENX_MTX):
            return "10x_mtx"
        # Cell Ranger often nests the matrix one level down
        for sub in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix", "outs"):
            p = os.path.join(path, sub)
            if os.path.isdir(p) and any(m in os.listdir(p) for m in _TENX_MTX):
                return "10x_mtx"
        return None
    lower = path.lower()
    if lower.endswith(".h5ad"):
        return "h5ad"
    if lower.endswith(".h5"):
        return "h5"
    return None


def resolve_matrix_dir(path):
    """If a 10x matrix lives in a known subdirectory, return that subdirectory."""
    if not os.path.isdir(path):
        return path
    names = set(os.listdir(path))
    if any(m in names for m in _TENX_MTX):
        return path
    for sub in ("filtered_feature_bc_matrix", "raw_feature_bc_matrix", "outs"):
        p = os.path.join(path, sub)
        if os.path.isdir(p) and any(m in os.listdir(p) for m in _TENX_MTX):
            return p
        # outs/filtered_feature_bc_matrix
        p2 = os.path.join(p, "filtered_feature_bc_matrix")
        if os.path.isdir(p2) and any(m in os.listdir(p2) for m in _TENX_MTX):
            return p2
    return path


def load_samplesheet(path, input_mode="matrix", check_files=True, base_dir=None):
    """Parse and validate. Returns (samples, extra_columns)."""
    if input_mode not in ("matrix", "fastq"):
        raise SampleSheetError(f"input_mode must be 'matrix' or 'fastq', got {input_mode!r}")
    if not os.path.exists(path):
        raise SampleSheetError(f"Sample sheet not found: {path}")

    delim = _sniff_delimiter(path)
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        if not reader.fieldnames:
            raise SampleSheetError(f"{path} has no header row.")
        mapping, unmapped = _map_columns(reader.fieldnames)
        canonical = set(mapping.values())
        if "sample" not in canonical:
            raise SampleSheetError(
                f"No sample-name column found in {path}.\n"
                f"Columns present: {reader.fieldnames}\n"
                f"Name one of them 'sample' (or sample_id / name / library).")
        need = "matrix" if input_mode == "matrix" else "fastq_1"
        if need not in canonical:
            hint = ("'matrix' (or path / counts / h5ad / directory)" if input_mode == "matrix"
                    else "'fastq_1' (or R1 / read1)")
            raise SampleSheetError(
                f"input_mode is {input_mode!r} but {path} has no {need!r} column.\n"
                f"Columns present: {reader.fieldnames}\nName one of them {hint}.")
        rows = list(reader)

    if not rows:
        raise SampleSheetError(f"{path} has a header but no data rows.")

    base_dir = base_dir or os.path.dirname(os.path.abspath(path))

    def resolve(p):
        p = (p or "").strip()
        if not p:
            return ""
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(base_dir, p))

    extra_cols = [c for c in mapping.values()
                  if c not in ("sample", "matrix", "fastq_1", "fastq_2")] + list(unmapped)

    samples, errors = OrderedDict(), []
    for i, row in enumerate(rows, start=2):
        rec = {}
        for raw, canon in mapping.items():
            rec[canon] = (row.get(raw) or "").strip()
        for raw in unmapped:
            rec[raw] = (row.get(raw) or "").strip()

        name = rec.get("sample", "")
        if not name:
            errors.append(f"line {i}: empty sample name")
            continue
        if not _SAFE_SAMPLE.match(name):
            errors.append(f"line {i}: sample name {name!r} must start alphanumeric and "
                          f"contain only letters, digits, dot, dash, underscore")
            continue
        if name in samples:
            errors.append(f"line {i}: duplicate sample name {name!r} "
                          f"(single-cell samples are not merged across rows)")
            continue

        entry = {"name": name}
        if input_mode == "matrix":
            m = resolve(rec.get("matrix", ""))
            if not m:
                errors.append(f"line {i} ({name}): no matrix path given")
                continue
            if check_files:
                if not os.path.exists(m):
                    errors.append(f"line {i} ({name}): matrix not found: {m}")
                    continue
                kind = detect_matrix_kind(m)
                if kind is None:
                    errors.append(
                        f"line {i} ({name}): could not identify {m} as a count matrix. "
                        f"Expected a directory containing matrix.mtx(.gz), or a .h5 / .h5ad file")
                    continue
                entry["kind"] = kind
                entry["matrix"] = resolve_matrix_dir(m) if kind == "10x_mtx" else m
            else:
                entry["kind"] = detect_matrix_kind(m) or "unknown"
                entry["matrix"] = m
        else:
            f1, f2 = resolve(rec.get("fastq_1", "")), resolve(rec.get("fastq_2", ""))
            if not f1:
                errors.append(f"line {i} ({name}): no FASTQ given")
                continue
            if not f2:
                errors.append(
                    f"line {i} ({name}): droplet single-cell data needs both reads "
                    f"(R1 carries the cell barcode + UMI, R2 the cDNA); fastq_2 is empty")
                continue
            if check_files:
                for f in (f1, f2):
                    if not os.path.exists(f):
                        errors.append(f"line {i} ({name}): FASTQ not found: {f}")
            entry["fastq_1"], entry["fastq_2"] = f1, f2

        for c in extra_cols:
            entry[c] = rec.get(c, "")
        samples[name] = entry

    if errors:
        raise SampleSheetError("Sample sheet problems in {}:\n  - {}".format(
            path, "\n  - ".join(errors)))

    # the same input assigned to two samples is nearly always a copy-paste error
    seen = defaultdict(list)
    for s in samples.values():
        for key in ("matrix", "fastq_1", "fastq_2"):
            if s.get(key):
                seen[s[key]].append(s["name"])
    dupes = {f: sorted(set(n)) for f, n in seen.items() if len(set(n)) > 1}
    if dupes:
        msg = "\n  - ".join(f"{f} used by {n}" for f, n in dupes.items())
        raise SampleSheetError(f"The same input is assigned to multiple samples:\n  - {msg}")

    return samples, extra_cols


def summarise(samples, extra_cols, input_mode="matrix"):
    lines = [f"{len(samples)} samples ({input_mode} input)"]
    if input_mode == "matrix":
        kinds = defaultdict(int)
        for s in samples.values():
            kinds[s.get("kind", "unknown")] += 1
        lines.append("  formats: " + ", ".join(f"{k}={v}" for k, v in sorted(kinds.items())))
    if extra_cols:
        lines.append(f"  metadata columns: {', '.join(extra_cols)}")
        for c in extra_cols:
            vals = sorted({s.get(c, "") for s in samples.values() if s.get(c, "")})
            if 0 < len(vals) <= 12:
                lines.append(f"    {c}: {', '.join(vals)}")
    return "\n".join(lines)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        print("usage: python scripts/samplesheet.py <samplesheet> [matrix|fastq]")
        return 2
    mode = argv[2] if len(argv) > 2 else "matrix"
    try:
        samples, extra = load_samplesheet(argv[1], input_mode=mode, check_files=True)
    except SampleSheetError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return 1
    print(summarise(samples, extra, mode))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
