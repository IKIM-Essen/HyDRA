rule fastqc_before_trim:
    input:
        "data/illumina/{strain}_{lane}_{read}.fastq"
    output:
        html="results/reports/before_trim/{strain}_{lane}_{read}_fastqc.html",
        zip="results/reports/before_trim/{strain}_{lane}_{read}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{strain}_{lane}_{read}.log"
    threads: 1
    wrapper:
        "v0.80.1/bio/fastqc"

use rule fastqc_before_trim as fastqc_trimmed with:
    input:
        "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq.gz"
    output:
        html= "results/reports/trimmed/{strain}_{read}_fastqc.html",
        zip= "results/reports/trimmed/{strain}_{read}_fastqc.zip"
    log:
        "logs/fastqc/{strain}_{read}_trimmed.log"