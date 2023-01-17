from pathlib import Path

'''## faster than old script and removes more
rule fastq_filter:
    input:
        r1 = get_ill_rawR1,
        r2 = get_ill_rawR2#("R2")
    output:
        r1_filt = "results/preprocess_ill/{strain}/{strain}_{lane}_R1_filtered.fastq.gz",
        r2_filt = "results/preprocess_ill/{strain}/{strain}_{lane}_R2_filtered.fastq.gz"
    params:
        # mean quality; also possible to filter by median quality
        quality = 26
    log:
        "logs/fastq_filter/{strain}_{lane}.log"
    conda:
        "../envs/fastq_filter.yaml"
    shell:
        "fastq-filter -q {params.quality} -o {output.r1_filt} -o {output.r2_filt} {input.r1} {input.r2} "
        "2> {log}"

rule cutadapt:
    input:
        ["results/preprocess_ill/{strain}/{strain}_{lane}_R1_filtered.fastq.gz", "results/preprocess_ill/{strain}/{strain}_{lane}_R2_filtered.fastq.gz"]
        #["results/preprocess_ill/{strain}/{strain}_{lane}_R1_dedup.fastq", "results/preprocess_ill/{strain}/{strain}_{lane}_R2_dedup.fastq"]
    output:
        fastq1="results/preprocess_ill/{strain}/{strain}_{lane}_R1_trimmed.fastq.gz",
        fastq2="results/preprocess_ill/{strain}/{strain}_{lane}_R2_trimmed.fastq.gz",
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
        l1 = "results/preprocess_ill/{strain}/{strain}_L001_{read}_trimmed.fastq.gz",
        l2 = "results/preprocess_ill/{strain}/{strain}_L002_{read}_trimmed.fastq.gz"
    output:
        "results/preprocess_ill/{strain}/{strain}_{read}.fastq.gz"
    #params:
        #fastq_out = "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq"
    shell:
        "cat {input.l1} {input.l2} > {output}"
        #{params.fastq_out} && ""gzip {params.fastq_out}"

## still some duplicates in R2..
rule deduplicate:
    input:
        r1 = "results/preprocess_ill/{strain}/{strain}_R1.fastq.gz",
        r2 = "results/preprocess_ill/{strain}/{strain}_R2.fastq.gz"
    output:
        r1_dedup = "results/preprocess_ill/{strain}/{strain}_R1_trimmed.fastq.gz",
        r2_dedup = "results/preprocess_ill/{strain}/{strain}_R2_trimmed.fastq.gz"
    log:
        "logs/deduplicate/{strain}.log"
    script:
        "../scripts/deduplication.py"
'''

### alternative way

rule fastp:
    input:
        r1 = get_ill_rawR1,
        r2 = get_ill_rawR2
    output:
        r1 = "results/preprocess_ill/{strain}/{strain}_{lane}_R1_trimmed.fastq.gz",
        r2 = "results/preprocess_ill/{strain}/{strain}_{lane}_R2_trimmed.fastq.gz",
        html = "results/preprocess_ill/{strain}/{strain}_{lane}_fastp.html",
        json = "results/preprocess_ill/{strain}/{strain}_{lane}_fastp.json"
    log:
        "logs/fastp/{strain}_{lane}_fastp.log"
    params:
        ## -a: adapter sequence to trim
        ## -q, -u, -e: quality filter
        ## -f: trim front bases
        ## --cut_tail -M 30: cut from tail, if mean quality in window of 4 < -M score
        ## -l: minimal read length
        ## -D: deduplicate
        ## -x: trim poly-x tail
        extra="-a CTGTCTCTTATACACATCT -q 28 -u 4 -e 30 -f 16 --cut_tail -M 30 -l 65 -x -D"
    conda:
        "../envs/fastp.yaml"
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-h {output.html} -j {output.json} {params.extra} 2> {log}"

rule combine_fastp_lanes:
    input:
        l1 = "results/preprocess_ill/{strain}/{strain}_L001_{read}_trimmed.fastq.gz",
        l2 = "results/preprocess_ill/{strain}/{strain}_L002_{read}_trimmed.fastq.gz"
    output:
        "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq.gz"
    shell:
        "cat {input.l1} {input.l2} > {output}"

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


## überflüssig?
rule extract_fastqc_summary:
    input:
        zip = "results/reports/trimmed/{strain}_{read}_fastqc.zip"
    output:
        "results/reports/trimmed/{strain}_{read}_fastqc_data.txt"
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
        """