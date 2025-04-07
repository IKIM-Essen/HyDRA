if get_has_short_reads():

    rule fastqc_ill:
        input:
            "results/{date}/trimmed/fastp/{sample}.{read}.fastq.gz",
        output:
            html=temp("results/{date}/qc/fastqc/{sample}/ill_{sample}_{read}.html"),
            zip=temp("results/{date}/qc/fastqc/{sample}/ill_{sample}_{read}_fastqc.zip"),
        params:
            extra="--quiet",
        log:
            "logs/{date}/illumina/fastqc/{sample}_{read}.log",
        threads: 10
        resources:
            mem_mb=1024,
        wrapper:
            "v3.10.2/bio/fastqc"


if get_has_long_reads():

    rule NanoPlot:
        input:
            local(rules.chopper.output.trim_filt),
        output:
            txt=local("results/{date}/qc/nanoplot/{sample}/{sample}_NanoStats.txt"),
            html=local(
                "results/{date}/qc/nanoplot/{sample}/{sample}_NanoPlot-report.html"
            ),
        params:
            outdir=lambda wildcards, output: Path(output.txt).parent,
            extra="--huge -f svg",
        log:
            "logs/{date}/ont/nanoplot/{sample}.log",
        threads: 20
        conda:
            "../envs/nanoplot.yaml"
        shell:
            "NanoPlot --fastq {input} --prefix {wildcards.sample}_ -t {threads} "
            "{params.extra} --outdir {params.outdir} > {log} 2>&1"


rule NanoPlotRaw:
    input:
        local(rules.copy_fastq_ont.output.fastq),
    output:
        txt="results/{date}/qc/nanoplot/raw/{sample}/{sample}_NanoStats.txt",
        html="results/{date}/qc/nanoplot/raw/{sample}/{sample}_NanoPlot-report.html",
    params:
        outdir=lambda wildcards, output: Path(output.txt).parent,
        extra="--huge -f svg",
    log:
        "logs/{date}/ont/nanoplot/{sample}_raw.log",
    threads: 20
    conda:
        "../envs/nanoplot.yaml"
    shell:
        "NanoPlot --fastq {input} --prefix {wildcards.sample}_ -t {threads} "
        "{params.extra} --outdir {params.outdir} > {log} 2>&1"


rule multiqc:
    input:
        local(get_multiqc_input),
        config=get_multiqc_config(),
    output:
        "results/{date}/report/multiqc.html",
        "results/{date}/report/multiqc_data.zip",
    params:
        extra="--zip-data-dir",
        use_input_files_only=True,
    threads: 20
    log:
        "logs/{date}/qc/multiqc.log",
    wrapper:
        "v3.10.2/bio/multiqc"


"""
rule qc_summary:
    input:
        qc_csv=rules.qc_summary.output.csv,
        asbl=expand(
            "results/{{date}}/report_prerequisites/assembly/{sample}_megahit.log",
            sample=get_samples(),
        ),
        mapped=expand(
            "results/{{date}}/report_prerequisites/assembly/{sample}_reads_mapped.txt",
            sample=get_samples(),
        ),
    output:
        csv="results/{date}/out/report/all/assembly_summary.csv",
        vis_csv=temp("results/{date}/out/report/all/assembly_summary_visual.csv"),
    log:
        "logs/{date}/report/qc_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/qc_summary.py
"""
