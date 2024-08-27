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
        threads: 20
        resources:
            mem_mb=1024,
        wrapper:
            "v3.10.2/bio/fastqc"


if get_has_long_reads():

    rule NanoPlot_prefilt:
        input:
            get_ont_fastq,
        output:
            txt="results/{date}/qc/prefilt_nanoplot/{sample}/{sample}_NanoStats.txt",
            html="results/{date}/qc/prefilt_nanoplot/{sample}/{sample}_NanoPlot-report.html",
        params:
            outdir=lambda wildcards, output: Path(output.txt).parent,
            extra="--huge -f svg",
        log:
            "logs/{date}/ont/nanoplot_prefilt/{sample}.log",
        threads: 20
        conda:
            "../envs/nanoplot.yaml"
        shell:
            "NanoPlot --fastq {input} --prefix {wildcards.sample}_ -t {threads} "
            "{params.extra} --outdir {params.outdir} > {log} 2>&1"

    rule NanoPlot:
        input:
            rules.chopper.output.trim_filt,
        output:
            txt="results/{date}/qc/nanoplot/{sample}/{sample}_NanoStats.txt",
            html="results/{date}/qc/nanoplot/{sample}/{sample}_NanoPlot-report.html",
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


rule multiqc_prefilt:
    input:
        expand(
            [
                "results/{{date}}/qc/prefilt_nanoplot/{sample}/{sample}_NanoStats.txt",
                "results/{{date}}/qc/prefilt_nanoplot/{sample}/{sample}_NanoPlot-report.html",
            ],
            sample=get_samples(),
        ),
        config=get_multiqc_config(),
    output:
        html="results/{date}/report_prefilt/multiqc.html",
        arc="results/{date}/report_prefilt/multiqc_data.zip",
    params:
        extra="--zip-data-dir",
        use_input_files_only=True,
    threads: 20
    log:
        "logs/{date}/qc/multiqc.log",
    wrapper:
        "v3.10.2/bio/multiqc"


rule multiqc:
    input:
        get_multiqc_input,
        config=get_multiqc_config(),
    output:
        html="results/{date}/report/multiqc.html",
        arc="results/{date}/report/multiqc_data.zip",
    params:
        extra="--zip-data-dir",
        use_input_files_only=True,
    threads: 20
    log:
        "logs/{date}/qc/multiqc.log",
    wrapper:
        "v3.10.2/bio/multiqc"


rule unzip_multiqc_report:
    input:
        prefilt=rules.multiqc_prefilt.output.arc,
        postfilt=rules.multiqc.output.arc,
    output:
        prefilt="results/{date}/report_prefilt/multiqc_data.json",
        postfilt="results/{date}/report/multiqc_data.json",
    params:
        prefilt=lambda wildcards, output: Path(output.prefilt).parent,
        postfilt=lambda wildcards, output: Path(output.postfilt).parent,
    log:
        "logs/{date}/qc/multiqc_unzip.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "unzip {input.prefilt} -d {params.prefilt} && "
        "unzip {input.postfilt} -d {params.postfilt} > {log} 2>&1"


rule qc_summary:
    input:
        prefilt=rules.unzip_multiqc_report.output.prefilt,
        postfilt=rules.unzip_multiqc_report.output.postfilt,
    output:
        csv_prefilt="results/{date}/report_prefilt/QC_all_prefilt.csv",
        csv_postfilt="results/{date}/report/QC_all_postfilt.csv",
        #vis_csv=temp("results/{date}/out/report/all/assembly_summary_visual.csv"),
    params:
        samples=get_samples(),
    log:
        "logs/{date}/report/qc_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/qc_summary.py"

