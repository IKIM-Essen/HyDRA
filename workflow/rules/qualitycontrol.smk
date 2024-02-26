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
        threads: 2
        resources:
            mem_mb=1024,
        wrapper:
            "v3.3.6/bio/fastqc"


'''
use rule fastqc_ill as fastqc_ont with:
    input:
        "results/{date}/trimmed/chopper/{sample}.fastq.gz",
    output:
        html=temp("results/{date}/qc/fastqc/{sample}/ont_{sample}.html"),
        zip=temp("results/{date}/qc/fastqc/{sample}/ont_{sample}_fastqc.zip"),
    log:
        "logs/{date}/ont/fastqc/{sample}.log",
'''

if get_has_long_reads():
    rule NanoPlot:
        input:
            rules.chopper.output.trim_filt,
        output:
            txt = "results/{date}/qc/nanoplot/{sample}/{sample}_NanoStats.txt",
            html = "results/{date}/qc/nanoplot/{sample}/{sample}_NanoPlot-report.html"
        params:
            outdir = lambda wildcards, output: Path(output.txt).parent,
            extra = "--huge -f svg",
        log:
            "logs/{date}/ont/nanoplot/{sample}.log"
        threads: 4
        conda:
            "../envs/nanoplot.yaml"
        shell:
            "NanoPlot --fastq {input} --prefix {wildcards.sample}_ -t {threads} "
            "{params.extra} --outdir {params.outdir} > {log} 2>&1"


rule multiqc:
    input:
        get_multiqc_input,
    output:
        "results/{date}/report/multiqc.html",
        "results/{date}/report/multiqc_data.zip",
    params:
        extra="--zip-data-dir",
        use_input_files_only=True,
    log:
        "logs/{date}/qc/multiqc.log",
    wrapper:
        "v3.3.6/bio/multiqc"
