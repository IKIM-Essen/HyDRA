from pathlib import Path
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
    log:
        "logs/fastp/{strain}_{read}_fastp.log"
    conda:
        "../envs/fastp.yaml"
    shell:
        "cat {input.l1} {input.l2} > {output} 2> {log}"

