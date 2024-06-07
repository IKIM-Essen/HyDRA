# copy fastq files to local
if get_has_short_reads():

    rule copy_fastq_illumina:
        input:
            fastqs=get_ill_fastqs,
        output:
            fastqs=[
                "{0}{{date}}/{{sample}}_R1.fastq.gz".format(get_data_path_ill()),
                "{0}{{date}}/{{sample}}_R2.fastq.gz".format(get_data_path_ill()),
            ],
        params:
            indir=lambda wildcards, input: Path(input.fastqs[0]).parent,
            infiles=[
                lambda wildcards, input: Path(input.fastqs[0]).name,
                lambda wildcards, input: Path(input.fastqs[1]).name,
            ],
            outdir=lambda wildcards, output: Path(output.fastqs[0]).parent,
        log:
            "logs/{date}/copy_data/{sample}_ill.log",
        threads: 64
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.outdir} && "
            "(cd {params.indir} && "
            "tar cpfz - {input.infiles}) | "
            "(cd {params.outdir} ; tar xpfz - )) > {log} 2>&1"


if get_has_long_reads():

    rule copy_fastq_ont:
        input:
            fastq=get_ont_fastq,
        output:
            fastq="{0}{{date}}/{{sample}}.fastq.gz".format(get_data_path_ont()),
        params:
            indir=lambda wildcards, input: Path(input.fastq).parent,
            infile=lambda wildcards, input: Path(input.fastq).name,
            outdir=lambda wildcards, output: Path(output.fastq).parent,
        log:
            "logs/{date}/copy_data/{sample}_ont.log",
        threads: 64
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.outdir} && "
            "(cd {params.indir} && "
            "tar cpfz - {input.infile}) | "
            "(cd {params.outdir} ; tar xpfz - )) > {log} 2>&1"
