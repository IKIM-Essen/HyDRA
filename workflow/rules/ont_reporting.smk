rule pycoqc:
    input:
        get_summaryfile_by_run
    output:
        "results/reports/run_QC/{run}_pycoqc.html"
    log:
        "logs/pycoqc/{run}.log"
    conda:
        "../envs/pycoqc.yaml"
    shell:
        "pycoQC -f {input} -o {output} 2> {log}"#" 2>&1"

rule nanoQC:
    input:
        get_ont_reads_by_stage
    output:
        html = "results/reports/{stage}/{strain}/nanoQC.html",
        log = "results/reports/{stage}/{strain}/NanoQC.log"
    params:
        outdir = lambda wildcards, output: Path(output.html).parent
    #log file is created by the tool & moved to the log folder
    conda:
        "../envs/nanoqc.yaml"
    shell:
        "mkdir -p {params.outdir} && "
        "nanoQC -o {params.outdir} {input}"

rule NanoPlot:
    input:
        get_ont_reads_by_stage
    output:
        html = "results/reports/{stage}/{strain}/NanoPlot/{strain}_NanoPlot-report.html",
        txt = "results/reports/{stage}/{strain}/NanoPlot/{strain}_NanoStats.txt"
    params:
        outdir = lambda wildcards, output: Path(output.txt).parent,
        prefix = "{strain}_"
    log:
        "logs/nanoplot/{strain}_{stage}.log"
    conda:
        "../envs/nanoplot.yaml"
    shell:
        "NanoPlot --fastq {input} --prefix {params.prefix} --verbose --outdir {params.outdir} 2> {log}"#" 2>&1"

rule rename_ont_reports:
    input:
        nanostats = "results/reports/{stage}/{strain}/NanoPlot/{strain}_NanoStats.txt",
        nanoplot = "results/reports/{stage}/{strain}/NanoPlot/{strain}_NanoPlot-report.html",
        nanoqc = "results/reports/{stage}/{strain}/nanoQC.html",
        nanoqc_log = "results/reports/{stage}/{strain}/NanoQC.log"

    output:
        nanostats = "results/reports/{stage}/{strain}/{strain}_{stage}_NanoStats.txt",
        nanoplot = "results/reports/{stage}/{strain}/{strain}_{stage}_NanoPlot-report.html",
        nanoqc = "results/reports/{stage}/{strain}/{strain}_{stage}_nanoQC.html",
        nanoqc_log = "logs/nanoqc/{strain}_{stage}.log"
    shell:
        "mv {input.nanostats} {output.nanostats} && "
        "mv {input.nanoplot} {output.nanoplot} && "
        "mv {input.nanoqc} {output.nanoqc} && "
        "mv {input.nanoqc_log} {output.nanoqc_log}"
