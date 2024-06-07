from pathlib import Path

ASSEMBLY_TYPE = get_assembly_type()

if ASSEMBLY_TYPE == "hybrid":

    rule unicycler_hybrid:
        input:
            # R1 and R2 short reads:
            paired=rules.fastp.output.trimmed,
            # Long reads:
            long=rules.chopper.output.trim_filt,
        output:
            assembly=temp("results/{date}/assembly/{sample}/assembly.fasta"),
        log:
            "logs/{date}/assembly/unicycler_hybrid/{sample}.log",
        params:
            extra=" --min_fasta_length 500 ",
        threads: 64
        wrapper:
            "v3.10.2/bio/unicycler"

elif ASSEMBLY_TYPE == "short":

    rule unicycler_short:
        input:
            # R1 and R2 short reads:
            paired=rules.fastp.output.trimmed,
        output:
            assembly=temp(get_assembly),
        log:
            "logs/{date}/assembly/unicycler_short/{sample}.log",
        params:
            extra="--min_fasta_length 300 --keep 0",
        threads: 64
        wrapper:
            "v3.10.2/bio/unicycler"

elif ASSEMBLY_TYPE == "long":

    rule unicycler_long:
        input:
            # Long reads:
            long=rules.chopper.output.trim_filt,
        output:
            assembly=temp(get_assembly),
        log:
            "logs/{date}/assembly/unicycler_long/{sample}.log",
        params:
            extra="--min_fasta_length 300 --keep 0",
        threads: 64
        wrapper:
            "v3.10.2/bio/unicycler"


rule assembly_gz:
    input:
        get_assembly,
    output:
        fa_gz="results/{date}/out/assembly/{sample}.fa.gz",
    log:
        "logs/{date}/assembly/gz_{sample}.log",
    threads: 20
    conda:
        "../envs/unix.yaml"
    shell:
        "pigz -c {input} > {output.fa_gz} 2> {log}"


if ASSEMBLY_TYPE in ["hybrid", "short"]:

    rule map_ill_to_assembly:
        input:
            contigs=get_assembly,
            fastqs=rules.fastp.output.trimmed,
        output:
            "results/{date}/report_prerequisites/assembly/{sample}_short_reads_mapped.txt",
        threads: 64
        log:
            "logs/{date}/assembly/{sample}_mapping_short_reads.log",
        conda:
            "../envs/minimap2.yaml"
        shell:
            "(minimap2 -ax sr -t {threads} {input.contigs} {input.fastqs} | "
            "samtools view -c -F 4 --threads {threads} -o {output}) > {log} 2>&1"


if ASSEMBLY_TYPE in ["hybrid", "long"]:

    # returns number of primary mapped reads
    rule map_ont_to_assembly:
        input:
            contigs=get_assembly,
            fastq=rules.chopper.output.trim_filt,
        output:
            aln="results/{date}/report_prerequisites/assembly/{sample}_long_reads_mapped.txt",
        threads: 64
        log:
            "logs/{date}/assembly/{sample}_mapping_long_reads.log",
        conda:
            "../envs/minimap2.yaml"
        shell:
            "(minimap2 -ax map-ont -t {threads} {input.contigs} {input.fastq} | "
            "samtools view -c -F 2308 --threads {threads} -o {output}) > {log} 2>&1"


"""
rule assembly_summary:
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
        csv="results/{date}/output/report/all/assembly_summary.csv",
        vis_csv=temp("results/{date}/output/report/all/assembly_summary_visual.csv"),
    log:
        "logs/{date}/report/assembly_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assembly_summary.py
"""
