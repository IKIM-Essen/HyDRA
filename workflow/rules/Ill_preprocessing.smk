from pathlib import Path

rule fastqc_before_trim:
    input:
        "data/illumina/{strain}_{lane}_{read}.fastq"
    output:
        html="results/reports/before_trim/all_reports/{strain}_{lane}_{read}_fastqc.html",
        zip="results/reports/before_trim/fastqc/{strain}/{strain}_{lane}_{read}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{strain}_{lane}_{read}.log"
    threads: 1
    wrapper:
        "v0.80.1/bio/fastqc"

rule minirmd:
    input:
        fastq1 = get_Illread1_by_lane,
        fastq2 = get_Illread2_by_lane
    output:
        "results/minirmd/{strain}_{lane}_1","results/minirmd/{strain}_{lane}_2"
    params:
        d = "2", #mismatch_number
        gen_outfile = "results/minirmd/{strain}_{lane}"
    log:
        "logs/minirmd/{strain}_{lane}.log"
    conda:
        "../envs/minirmd.yaml"
    shell:
        "minirmd -i {input.fastq1} -f {input.fastq2} -o {params.gen_outfile} -d {params.d} > {log} 2>&1"

rule minirmd_renaming:
    input:
        fastq1 = "results/minirmd/{strain}_{lane}_1",
        fastq2 = "results/minirmd/{strain}_{lane}_2"
    output:
        fastq1="results/minirmd/{strain}_{lane}_R1_uniq.fastq",
        fastq2="results/minirmd/{strain}_{lane}_R2_uniq.fastq"
    shell:
        "mv {input.fastq1} {output.fastq1} && "
        "mv {input.fastq2} {output.fastq2}"

"""
rule cutadapt:
    input:
        ["results/minirmd/{strain}_{lane}_R1_uniq.fastq", "results/minirmd/{strain}_{lane}_R2_uniq.fastq"]
    output:
        fastq1="results/cutadapt/{strain}/{strain}_{lane}_R1_trimmed.fastq",
        fastq2="results/cutadapt/{strain}/{strain}_{lane}_R2_trimmed.fastq",
        qc="results/cutadapt/{strain}/{strain}_{lane}.qc.txt"
    params:
        adapters= "-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT", #lambda adapt: -a \"G{20}\" -A \"G{20}\"
        #maybe quality of 27 is better (cause of quality trimming algorithm)
        extra="--nextseq-trim=26 -q 26,26 -u 15 -U 15 --minimum-length 136" # -O 10 -u 12 -U 12 -q 20,20 
    log:
        "logs/cutadapt/{strain}_{lane}.log"
    threads: 4
    wrapper:
        "v0.80.1/bio/cutadapt/pe"

rule cutadapt_impro:
    input:
        ["old_filter_results/cutadapt/{strain}/{strain}_{lane}_R1_trimmed.fastq", "old_filter_results/cutadapt/{strain}/{strain}_{lane}_R2_trimmed.fastq"]
    output:
        fastq1="results/cutadapt/{strain}/{strain}_{lane}_R1_trimmed.fastq",
        fastq2="results/cutadapt/{strain}/{strain}_{lane}_R2_trimmed.fastq",
        qc="results/cutadapt/{strain}/{strain}_{lane}.qc.txt"
    params:
        adapters= "-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT", #lambda adapt: -a \"G{20}\" -A \"G{20}\"
        extra="--minimum-length 136" # -O 10 -u 12 -U 12 -q 20,20
    log:
        "logs/cutadapt/{strain}_{lane}.log"
    threads: 4
    wrapper:
        "v0.80.1/bio/cutadapt/pe"

"""

rule qfilter:
    input:
        fastq1 = get_Illread1_by_lane,
        fastq2 = get_Illread2_by_lane
    output:
        fastq1_filt = "results/preprocess_ill/{strain}/{strain}_{lane}_R1_filtered.fastq",
        fastq2_filt = "results/preprocess_ill/{strain}/{strain}_{lane}_R2_filtered.fastq",
        unpaired = "results/preprocess_ill/{strain}/{strain}_{lane}_unpaired.fastq"
    params:
        q = "26",
        p = "90" 
    log:
        "logs/qfilter/{strain}_{lane}.log"
    conda:
        "../envs/fastx_toolkit.yaml"
    script:
        "../scripts/qfiltering.py"

rule deduplicate:
    input:
        fastq1 = "results/preprocess_ill/{strain}/{strain}_{lane}_R1_filtered.fastq",
        fastq2 = "results/preprocess_ill/{strain}/{strain}_{lane}_R2_filtered.fastq"
    output:
        fastq1_dedup = "results/preprocess_ill/{strain}/{strain}_{lane}_R1_dedup.fastq",
        fastq2_dedup = "results/preprocess_ill/{strain}/{strain}_{lane}_R2_dedup.fastq"
    log:
        "logs/deduplicate/{strain}_{lane}.log"
    script:
        "../scripts/deduplication.py"

rule cutadapt:
    input:
        ["results/preprocess_ill/{strain}/{strain}_{lane}_R1_dedup.fastq", "results/preprocess_ill/{strain}/{strain}_{lane}_R2_dedup.fastq"]
    output:
        fastq1="results/preprocess_ill/{strain}/{strain}_{lane}_R1_trimmed.fastq",
        fastq2="results/preprocess_ill/{strain}/{strain}_{lane}_R2_trimmed.fastq",
        qc="results/preprocess_ill/{strain}/{strain}_{lane}.qc.txt"
    params:
        adapters= "-a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT", #lambda adapt: -a \"G{20}\" -A \"G{20}\"
        #maybe quality of 27 is better (cause of quality trimming algorithm)
        extra="--nextseq-trim=27 -q 27,27 -u 15 -U 15 --minimum-length 65"  #27 -q 27,27
    log:
        "logs/cutadapt/{strain}_{lane}.log"
    threads: 4
    wrapper:
        "v0.80.1/bio/cutadapt/pe"

rule combine_lanes:
    input:
        l1 = "results/preprocess_ill/{strain}/{strain}_L001_{read}_trimmed.fastq",
        l2 = "results/preprocess_ill/{strain}/{strain}_L002_{read}_trimmed.fastq"
    output:
        "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq.gz"
    params:
        fastq_out = "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq"
    shell:
        "cat {input.l1} {input.l2} > {params.fastq_out} && "
        "gzip {params.fastq_out}"

rule fastqc_trimmed:
    input:
        "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq.gz"
    output:
        html="results/reports/trimmed/all_reports/{strain}_{read}_trimmed_fastqc.html",
        zip="results/reports/trimmed/fastqc/{strain}/{strain}_{read}_trimmed_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{strain}_{read}_trimmed.log"
    threads: 1
    wrapper:
        "v0.80.1/bio/fastqc"

rule extract_fastqc_summary:
    input:
        zip = "results/reports/trimmed/fastqc/{strain}/{strain}_{read}_trimmed_fastqc.zip"
    output:
        "results/reports/trimmed/fastqc/{strain}/{strain}_{read}_fastqc_data.txt"
    params:
        outdir = lambda wildcards, input: Path(input.zip).parent,
        zipdir = "results/reports/trimmed/fastqc/{strain}/{strain}_{read}_trimmed_fastqc/"
    shell:
        "unzip -d {params.outdir} {input.zip} && "
        "mv {params.zipdir}/fastqc_data.txt {output} && "
        "rm -r {params.zipdir}"

rule coverage_ill:
    input:
        "results/reports/trimmed/fastqc/{strain}/{strain}_R1_fastqc_data.txt"
    output:
        "results/reports/trimmed/fastqc/{strain}_coverage.txt"
    params:
        tech = "illumina",
        genomesize = get_genome_size
    conda:
        "../envs/watstats.yaml"
    script:
        "../scripts/coverage.py"

rule coverages_ill:
    input:
        expand("results/reports/trimmed/fastqc/{strain}_coverage.txt", strain=get_all_strain_ids())
    output:
        outfile = "results/reports/trimmed/all_reports/ill_coverages.txt"
    params:
        indir = "results/reports/trimmed/fastqc/" #lambda wildcards, output: Path(output.outfile).parent
    shell:
        "cat {input} > {output} && "
        "rm {params.indir}*coverage.txt"