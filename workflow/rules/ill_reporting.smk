rule fastqc_before_trim:
    input:
        "data/illumina/{strain}_{lane}_{read}.fastq"
        #get_ill_reads_by_stage
    output:
        html="results/reports/before_trim/{strain}_{lane}_{read}_fastqc.html",
        zip="results/reports/before_trim/{strain}_{lane}_{read}_fastqc.zip"
        #html="results/reports/{stage}/all_reports/{strain}_{lane}_{read}_fastqc.html",
        #zip="results/reports/{stage}/fastqc/{strain}/{strain}_{lane}_{read}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{strain}_{lane}_{read}.log"
    threads: 1
    wrapper:
        "v0.80.1/bio/fastqc"

rule fastqc_trimmed:
    input:
        "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq.gz"
        #get_ill_reads_by_stage
    output:
        html= "results/reports/trimmed/{strain}_{read}_fastqc.html",
        zip= "results/reports/trimmed/{strain}_{read}_fastqc.zip"
    params: "--quiet"
    log:
        "logs/fastqc/{strain}_{read}_trimmed.log"
    threads: 1
    wrapper:
        "v0.80.1/bio/fastqc"