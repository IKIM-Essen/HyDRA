# RAW_DATA_PATH = get_data_path_ill()


# copy fastq files to local
if get_has_short_reads():

    rule copy_fastq_illumina:
        input:
            fastqs=get_ill_fastqs,
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
            "cp -v {input.fastqs[0]} {output.fastqs[0]} && "
            "cp -v {input.fastqs[1]} {output.fastqs[1]}) > {log} 2>&1"


if get_has_long_reads():

    rule copy_fastq_ont:
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
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.outdir} && "
            "cp -v {input} {output.fastq}) > {log} 2>&1"
