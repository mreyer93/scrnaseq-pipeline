"""The analysis stages. Each writes an .h5ad so a long run can resume part-way and each
step's output can be inspected independently."""

_common = dict(
    figures_dir = FIGURES_DIR,
    tables_dir = TABLES_DIR,
)


rule qc:
    """Load every sample, compute QC metrics, filter cells and genes."""
    input:
        matrices = [SAMPLES[s]["matrix"] for s in SAMPLE_NAMES] if INPUT_MODE == "matrix"
                   else expand(join(OUTDIR, "quant", "{sample}", "counts.h5ad"),
                               sample=SAMPLE_NAMES),
    output:
        h5ad = join(OUTDIR, "h5ad", "01_qc.h5ad"),
        metrics = metrics_path("qc"),
    log: join(OUTDIR, "logs", "qc.log")
    params:
        samples = lambda w: _samples_for_stage(),
        extra_cols = EXTRA_COLS,
        qc = config["qc"],
        metrics_file = metrics_path("qc"),
        **_common
    conda: "../../envs/python.yml"
    script: "../../scripts/qc.py"


def _samples_for_stage():
    """Sample dict pointing at whichever matrices this run actually uses."""
    if INPUT_MODE == "matrix":
        return {k: dict(v) for k, v in SAMPLES.items()}
    out = {}
    for s in SAMPLE_NAMES:
        d = dict(SAMPLES[s])
        d["matrix"] = join(OUTDIR, "quant", s, "counts.h5ad")
        d["kind"] = "h5ad"
        out[s] = d
    return out


rule doublets:
    """Per-sample Scrublet doublet detection."""
    input:  h5ad = join(OUTDIR, "h5ad", "01_qc.h5ad")
    output:
        h5ad = join(OUTDIR, "h5ad", "02_doublets.h5ad"),
        metrics = metrics_path("doublets"),
    log: join(OUTDIR, "logs", "doublets.log")
    params:
        doublets = config["doublets"],
        metrics_file = metrics_path("doublets"),
        **_common
    conda: "../../envs/python.yml"
    script: "../../scripts/doublets.py"


rule cluster:
    """Normalise, select features, reduce dimensions, integrate, cluster, embed."""
    input:  h5ad = join(OUTDIR, "h5ad", "02_doublets.h5ad")
    output:
        h5ad = join(OUTDIR, "h5ad", "03_clustered.h5ad"),
        metrics = metrics_path("cluster"),
    log: join(OUTDIR, "logs", "cluster.log")
    params:
        cluster = config["cluster"],
        integration = config["integration"],
        metrics_file = metrics_path("cluster"),
        **_common
    conda: "../../envs/python.yml"
    script: "../../scripts/cluster.py"


rule markers:
    """Marker genes per cluster and optional cell-type scoring."""
    input:  h5ad = join(OUTDIR, "h5ad", "03_clustered.h5ad")
    output:
        h5ad = join(OUTDIR, "h5ad", "04_annotated.h5ad"),
        markers = join(TABLES_DIR, "markers.tsv"),
        metrics = metrics_path("markers"),
    log: join(OUTDIR, "logs", "markers.log")
    params:
        markers = config["markers"],
        annotation = config["annotation"],
        metrics_file = metrics_path("markers"),
        **_common
    conda: "../../envs/python.yml"
    script: "../../scripts/markers.py"


_report_params = dict(
    scripts_dir = SCRIPTS_DIR,
    outdir = OUTDIR,
    project_name = config["project_name"],
    figures_dir = FIGURES_DIR,
    tables_dir = TABLES_DIR,
    metrics_dir = METRICS_DIR,
)


rule report_html:
    input:
        h5ad = join(OUTDIR, "h5ad", "04_annotated.h5ad"),
        markers = join(TABLES_DIR, "markers.tsv"),
        metrics = metrics_path("markers"),
    output: join(OUTDIR, "report", "scrnaseq_report.html")
    log: join(OUTDIR, "logs", "report_html.log")
    params: rmd = join(SCRIPTS_DIR, "scrnaseq_report.Rmd"), format = "html_document",
            **_report_params
    conda: "../../envs/r.yml"
    script: "../../scripts/render_report.R"


rule report_pdf:
    input:
        h5ad = join(OUTDIR, "h5ad", "04_annotated.h5ad"),
        markers = join(TABLES_DIR, "markers.tsv"),
        metrics = metrics_path("markers"),
    output: join(OUTDIR, "report", "scrnaseq_report.pdf")
    log: join(OUTDIR, "logs", "report_pdf.log")
    params: rmd = join(SCRIPTS_DIR, "scrnaseq_report.Rmd"), format = "pdf_document",
            **_report_params
    conda: "../../envs/r.yml"
    script: "../../scripts/render_report.R"
