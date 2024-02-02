# RAW_DATA_PATH = get_data_path_ill()


# copy fastq files to local
rule copy_fastq_illumina:
    input:
        get_ill_fastqs,
    output:
        fastqs=[
            "{0}{{date}}/{{sample}}_R1.fastq.gz".format(
                get_data_path_ill()
            ),
            "{0}{{date}}/{{sample}}_R2.fastq.gz".format(
                get_data_path_ill()
            ),
        ],
    params:
        outdir=lambda wildcards, output: Path(output.fastqs[0]).parent,
    log:
        "logs/{date}/copy_data/{sample}_ill.log",
    conda:
        "../envs/unix.yaml"
    shell:
        "(mkdir -p {params.outdir} && "
        "cp -v -t {params.outdir} {input}) > {log} 2>&1"


use rule copy_fastq_illumina as copy_fastq_ont with:
    input:
        get_ont_fastq,
    output:
        fastq="{0}{{date}}/{{sample}}.fastq.gz".format(
            get_data_path_ont()
        ),
    params:
        outdir=lambda wildcards, output: Path(output.fastq).parent,
    log:
        "logs/{date}/copy_data/{sample}_ont.log",
