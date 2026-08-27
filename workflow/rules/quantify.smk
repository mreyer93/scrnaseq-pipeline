"""FASTQ -> count matrix with simpleaf (piscem index + alevin-fry quantification).

This is the same default nf-core/scrnaseq uses. Cell Ranger is deliberately not the
default here: it is 10x's proprietary tool under their EULA, cannot be redistributed in
a conda environment, and needs a manual download and licence acceptance. simpleaf is
open, fast, and benchmarks comparably for gene-level quantification.

Only included when input_mode is 'fastq'. In 'matrix' mode (the local path) the pipeline
starts from existing count matrices and none of this runs.
"""


rule simpleaf_index:
    """Build a splici (spliced + intronic) index once; reused by every sample."""
    input:
        genome = config.get("reference", {}).get("genome_fasta", ""),
        gtf = config.get("reference", {}).get("gtf", ""),
    output: directory(join(OUTDIR, "quant", "_index"))
    log: join(OUTDIR, "logs", "simpleaf_index.log")
    threads: config["threads"]["quant"]
    params:
        rlen = config.get("reference", {}).get("read_length", 91),
    conda: "../../envs/environment.yml"
    shell:
        """
        export ALEVIN_FRY_HOME="${{ALEVIN_FRY_HOME:-$PWD/.alevin_fry_home}}"
        mkdir -p "$ALEVIN_FRY_HOME"
        simpleaf set-paths >> {log} 2>&1 || true
        simpleaf index \
            --output {output} \
            --fasta {input.genome} \
            --gtf {input.gtf} \
            --rlen {params.rlen} \
            --threads {threads} >> {log} 2>&1
        """


rule simpleaf_quant:
    input:
        index = join(OUTDIR, "quant", "_index"),
        r1 = lambda w: SAMPLES[w.sample]["fastq_1"],
        r2 = lambda w: SAMPLES[w.sample]["fastq_2"],
    output:
        outdir = directory(join(OUTDIR, "quant", "{sample}", "af_quant")),
        h5ad = join(OUTDIR, "quant", "{sample}", "counts.h5ad"),
    log: join(OUTDIR, "logs", "simpleaf_quant", "{sample}.log")
    threads: config["threads"]["quant"]
    params:
        chemistry = config.get("chemistry", "10xv3"),
        resolution = config.get("resolution_strategy", "cr-like"),
        workdir = join(OUTDIR, "quant", "{sample}"),
    conda: "../../envs/environment.yml"
    shell:
        """
        export ALEVIN_FRY_HOME="${{ALEVIN_FRY_HOME:-$PWD/.alevin_fry_home}}"
        mkdir -p "$ALEVIN_FRY_HOME"
        simpleaf quant \
            --reads1 {input.r1} \
            --reads2 {input.r2} \
            --index {input.index}/index \
            --chemistry {params.chemistry} \
            --resolution {params.resolution} \
            --expected-ori fw \
            --t2g-map {input.index}/index/t2g_3col.tsv \
            --unfiltered-pl \
            --threads {threads} \
            --output {params.workdir} > {log} 2>&1

        python3 {workflow.basedir}/../scripts/alevin_to_h5ad.py \
            {output.outdir} {output.h5ad} >> {log} 2>&1
        """
