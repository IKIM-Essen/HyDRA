# fastp in paired-end mode for Illumina paired-end data
if get_has_short_reads():

    rule fastp:
        input:
            sample=rules.copy_fastq_illumina.output.fastqs,
        output:
            trimmed=temp(
                [
                    "results/{date}/trimmed/fastp/{sample}.1.fastq.gz",
                    "results/{date}/trimmed/fastp/{sample}.2.fastq.gz",
                ]
            ),
            html="results/{date}/trimmed/fastp/{sample}.html",
            json="results/{date}/trimmed/fastp/{sample}.fastp.json",
        params:
            adapters=get_adapters(),
            extra="--overrepresentation_analysis -e {phred} -l {minlen} -u 20 -q 20".format(
                phred=(config["quality_criteria"]["illumina"]["min_quality"]),
                minlen=(config["quality_criteria"]["illumina"]["min_length"]),
            ),
        log:
            "logs/{date}/illumina/fastp/{sample}.log",
        threads: 32
        wrapper:
            "v3.10.2/bio/fastp"


if get_has_long_reads():

    rule porechop_abi:
        input:
            rules.copy_fastq_ont.output.fastq,
        output:
            adapt_trim=temp("results/{date}/trimmed/porechop/{sample}.fastq.gz"),
        log:
            "logs/{date}/ont/porechop/{sample}.log",  # chnage to qc/ont/
        threads: 32
        conda:
            "../envs/porechop.yaml"
        shell:
            "porechop_abi -abi -i {input} --threads {threads} "
            "-o {output.adapt_trim} > {log} 2>&1"

    # Prowler testen
    # faster version of NanoFilt
    rule chopper:
        input:
            rules.porechop_abi.output.adapt_trim,
        output:
            trim_filt=temp("results/{date}/trimmed/chopper/{sample}.fastq.gz"),
        params:
            quality=config["quality_criteria"]["ont"]["min_quality"],
            minlen=config["quality_criteria"]["ont"]["min_length"],
        log:
            "logs/{date}/ont/chopper/{sample}.log",
        threads: 32
        conda:
            "../envs/chopper.yaml"
        shell:
            "(gunzip -c {input} | "
            "chopper -q {params.quality} -l {params.minlen} --threads {threads} | "
            "gzip > {output}) > {log} 2>&1"
