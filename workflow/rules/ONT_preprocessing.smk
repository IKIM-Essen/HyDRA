rule pycoqc:
    input:
        get_summaryfile_by_run
    output:
        "results/reports/pycoqc/{run}.html"
    log:
        "logs/pycoqc/{run}.log"
    conda:
        "../envs/pycoqc.yaml"
    shell:
        "pycoQC -f {input} -o {output} > {log} 2>&1"
   
rule combine_nanopore_fastqs:
    input:
        get_nanopore_reads
    output:
        "data/ont/{strain}.fastq.gz"
    log:
        "logs/com_ont_fastqs/{strain}.log"
    script:
        "../scripts/combine_ont_fastqs.py"

rule nanoQC:
    input:
        get_reads_by_stage
    output:
        "results/reports/{stage}/nanoqc/{strain}/nanoQC.html"
    params:
        "results/reports/{stage}/nanoqc/{strain}/"
    log:
        "logs/nanoqc/{strain}_{stage}.log"
    conda:
        "../envs/nanoqc.yaml"
    shell:
        "mkdir -p {params} && "
        "nanoQC -o {params} {input} > {log} 2>&1"

rule NanoPlot:
    input:
        get_reads_by_stage
    output:
        "results/reports/{stage}/nanoplot/{strain}/NanoPlot-report.html",
        "results/reports/{stage}/nanoplot/{strain}/NanoStats.txt"
    params:
        "results/reports/{stage}/nanoplot/{strain}/"
    log:
        "logs/nanoplot/{strain}_{stage}.log"
    conda:
        "../envs/nanoplot.yaml"
    shell:
        "NanoPlot --fastq {input} --outdir {params} > {log} 2>&1"

rule rename_reports:
    input:
        nanoplot = "results/reports/{stage}/nanoplot/{strain}/NanoPlot-report.html",
        nanoqc = "results/reports/{stage}/nanoqc/{strain}/nanoQC.html"
    output:
        nanoplot = "results/reports/{stage}/all_reports/{strain}_{stage}_NanoPlot-report.html",
        nanoqc = "results/reports/{stage}/all_reports/{strain}_{stage}_nanoQC.html"
    shell:
        "cp {input.nanoplot} {output.nanoplot} && "
        "cp {input.nanoqc} {output.nanoqc}"

rule nanofilt_lite:
    input:
        "data/ont/{strain}.fastq.gz"
    params:
        length = "5000",
        quality = "11"
    output:
        "results/nanofilt/{stage}/{strain}_flite.fastq.gz"
    log:
        "logs/nanofilt/{strain}_{stage}.log"
    conda:
        "../envs/nanofilt.yaml"
    shell:
        "gunzip -c {input} | "
        "NanoFilt -q {params.quality} -l {params.length} | "
        "gzip > {output} 2> {log}"

rule porechop:
    input:
        "results/nanofilt/before_trim/{strain}_flite.fastq.gz"
    output:
        #must be a file not a folder
        "results/porechop/{strain}_{stage}.fastq"
    log:
        "logs/porechop/{strain}_{stage}.log"
    conda:
        "../envs/porechop.yaml"
    shell:
        "porechop -i {input} -o {output} > {log} 2>&1"

rule gzip:
    input:
        "results/porechop/{strain}_{stage}.fastq"
    output:
        "results/porechop/{strain}_{stage}.fastq.gz"
    shell:
        "gzip {input}"

rule nanofilt:
    input:
        "results/porechop/{strain}_{stage}.fastq.gz"
    params:
        length = "8000", #8000
        quality = get_ont_quality_filter
    output:
        "results/nanofilt/{stage}/{strain}_{stage}.fastq.gz"
    log:
        "logs/nanofilt/{strain}_{stage}.log"
    conda:
        "../envs/nanofilt.yaml"
    shell:
        "gunzip -c {input} | "
        "NanoFilt -q {params.quality} -l {params.length} | "
        "gzip > {output} 2> {log}"

rule coverage_ont:
    input:
        "results/reports/trimmed/nanoplot/{strain}/NanoStats.txt"
    output:
        "results/reports/trimmed/nanoplot/{strain}_coverage.txt"
    params:
        tech = "ont",
        genomesize = get_genome_size
    conda:
        "../envs/watstats.yaml"
    script:
        "../scripts/coverage.py"

rule coverages_ont:
    input:
        expand("results/reports/trimmed/nanoplot/{strain}_coverage.txt", strain=get_all_strain_ids())
    output:
        outfile = "results/reports/trimmed/all_reports/ont_coverages.txt"
    params:
        indir = "results/reports/trimmed/nanoplot/" #lambda wildcards, output: Path(output.outfile).parent
    shell:
        "cat {input} > {output} && "
        "rm {params.indir}*coverage.txt"