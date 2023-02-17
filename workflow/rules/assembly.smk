from pathlib import Path

rule unicycler:
    input:
        # R1 and R2 short reads:
        paired_1 = "results/preprocess_ill/{strain}/{strain}_R1_trimmed.fastq.gz", #"results/cutadapt/{strain}/{strain}_R1_trimmed.fastq.gz",
        paired_2 = "results/preprocess_ill/{strain}/{strain}_R2_trimmed.fastq.gz", #"results/cutadapt/{strain}/{strain}_R2_trimmed.fastq.gz",
        # Long reads:
        long = "results/nanofilt/trimmed/{strain}_trimmed.fastq.gz"
    output:
        fasta = "results/assembly/{strain}/{strain}_assembly.fasta",
        gfa = "results/assembly/{strain}/{strain}_assembly.gfa"
    log:
        "logs/unicycler/{strain}.log" #automated log file - we need to mv it
    params:
        extra = "--min_fasta_length 500",
        outdir = "results/assembly/{strain}/unicycler/",
        log_old = "results/assembly/{strain}/unicycler/unicycler.log", # path.join(outdir, "unicycler.log")
        fasta_old = "results/assembly/{strain}/unicycler/assembly.fasta",
        gfa_old = "results/assembly/{strain}/unicycler/assembly.gfa"
    conda:
        "../envs/unicycler.yaml"
    shell:
        "unicycler -1 {input.paired_1} -2 {input.paired_2} -l {input.long} {params.extra} -o {params.outdir} && "#2> {log}"
        "mv {params.log_old} {log} && "
        "mv {params.fasta_old} {output.fasta} && "
        "mv {params.gfa_old} {output.gfa}"

"""
# use part of spades for unicycler
rule spades:
    input:
        r1 = "results/cutadapt/{strain}/{strain}_R1_trimmed.fastq.gz",
        r2 = "results/cutadapt/{strain}/{strain}_R2_trimmed.fastq.gz",
        ont = "results/nanofilt/trimmed/{strain}_trimmed.fastq.gz"
    output:
        scaf = "results/spades/{strain}/scaffolds.fasta"
    log:
        "logs/spades/{strain}.log"
    params:
        isolate = "--isolate",
        outdir = lambda wildcards, output: Path(output.scaf).parent
    conda:
        "../envs/spades.yaml"
    shell:
        "spades.py --pe1-1 {input.r1} --pe1-2 {input.r2} {params.isolate} "
        "--nanopore {input.ont} -o {params.outdir}/ > {log} 2>&1"

rule rename_scaffolds:
    input:
        "results/spades/{strain}/scaffolds.fasta"
    output:
        "results/spades/assembly/{strain}/{strain}_scaffolds.fasta"
    shell:
        "cp {input} {output}"

rule rename_assembly:
    input:
        "results/unicycler/{strain}/assembly.fasta"
    output:
        "results/assembly/{strain}/{strain}_assembly.fasta"
    shell:
        "cp {input} {output}"

rule canu:
    input:
        "results/nanofilt/trimmed/{strain}_trimmed.fastq.gz"
    output:
        outfile = "results/canu/{strain}/{strain}.contigs.fasta"
    params:
        size = get_genome_size,
        outdir = lambda wildcards, output: Path(output.outfile).parent
    log:
        "logs/canu/{strain}.log"
    threads: 16
    conda:
        "../envs/canu.yaml"
    shell:
        "canu -p {wildcards.strain} -d {params.outdir}/ genomeSize={params.size}k -nanopore {input} > {log} 2>&1"

rule bwa_index:
    input:
        "results/canu/{strain}/{strain}.contigs.fasta"
    output:
        "results/canu/{strain}/{strain}.contigs.fasta.amb",
        "results/canu/{strain}/{strain}.contigs.fasta.ann",
        "results/canu/{strain}/{strain}.contigs.fasta.bwt",
        "results/canu/{strain}/{strain}.contigs.fasta.pac",
        "results/canu/{strain}/{strain}.contigs.fasta.sa"
    log:
        "logs/bwa_index/{strain}.log"
    params:
        prefix="results/canu/{strain}/{strain}.contigs.fasta",
        algorithm="bwtsw"
    wrapper:
        "0.84.0/bio/bwa/index"

rule bwa_mem:
    input:
        reads=["results/cutadapt/{strain}/{strain}_R1_trimmed.fastq.gz", "results/cutadapt/{strain}/{strain}_R2_trimmed.fastq.gz"],
        ind = "results/canu/{strain}/{strain}.contigs.fasta.amb"
    output:
        "results/bwa/{strain}/{strain}.bam"
    log:
        "logs/bwa_mem/{strain}.log"
    params:
        index="results/canu/{strain}/{strain}.contigs.fasta",
        extra="",
        sorting="none",  # Can be 'none', 'samtools' or 'picard'.
        sort_order="queryname",  # Can be 'queryname' or 'coordinate'.
        sort_extra="",  # Extra args for samtools/picard.
        tmp_dir=""  # Path to temp dir. (optional)
    threads: 8
    wrapper:
        "0.84.0/bio/bwa/mem"

rule samtools_sort:
    input:
        "results/bwa/{strain}/{strain}.bam"
    output:
        "results/bwa/{strain}/{strain}.sorted.bam"
    params:
        extra = "-m 4G",
        tmp_dir = "/tmp/"
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "0.84.0/bio/samtools/sort"

rule samtools_index:
    input:
        "results/bwa/{strain}/{strain}.sorted.bam"
    output:
        "results/bwa/{strain}/{strain}.sorted.bam.bai"
    log:
        "logs/samtools_index/{strain}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "0.84.0/bio/samtools/index"

rule pilon:
    input:
        genome = "results/canu/{strain}/{strain}.contigs.fasta",
        frags = "results/bwa/{strain}/{strain}.sorted.bam",
        ind_frags = "results/bwa/{strain}/{strain}.sorted.bam.bai"
    output:
        outfile = "results/pilon/{strain}/{strain}.fasta"
    params:
        outdir = lambda wildcards, output: Path(output.outfile).parent
    log:
        "logs/pilon/{strain}.log"
    conda:
        "../envs/pilon.yaml"
    shell:
        "pilon --genome {input.genome} --frags {input.frags} --output {wildcards.strain} "
        "--outdir {params.outdir}/ > {log} 2>&1"

rule flye:
    input:
        "results/nanofilt/trimmed/{strain}_trimmed.fastq.gz"
    output:
        outfile = "results/flye/{strain}/assembly.fasta"
    params:
        outdir = lambda wildcards, output: Path(output.outfile).parent
    shell:
        "flye --nano-raw {input} --genome-size 4.6m -o {params.outdir}/"

rule medaka:
    input:
        "results/nanofilt/trimmed/{strain}_trimmed.fastq.gz"
    output:
        "results/medaka/{strain}_list_models.txt"
    conda:
        "../envs/medaka.yaml"
    shell:
        "medaka tools list_models > {output}"
"""