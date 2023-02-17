
'''rule combine_ont_fastqs:
    input:
        get_nanopore_reads
    output:
        "data/ont/{strain}.fastq.gz"
    log:
        "logs/com_ont_fastqs/{strain}.log"
    script:
        "../scripts/ont_combine_rawData.py"'''

rule nanofilt_lite:
    input:
        get_nanopore_reads
    params:
        length = "5000",
        quality = "11"
    output:
        "results/nanofilt/{stage}/{strain}_flite.fastq.gz"
    log:
        "logs/nanofilt/{strain}_{stage}.log"
    conda:
        "../envs/nanofilt.yaml"
    shell:
        "gunzip -c {input} | "
        "NanoFilt -q {params.quality} -l {params.length} --logfile {log} | "
        "gzip > {output}"

rule porechop:
    input:
        "results/nanofilt/before_trim/{strain}_flite.fastq.gz"
    output:
        #must be a file not a folder
        "results/porechop/{strain}_{stage}.fastq.gz"
    log:
        "logs/porechop/{strain}_{stage}.log"
    conda:
        "../envs/porechop.yaml"
    shell:
        "porechop -i {input} -o {output} 2> {log}"

use rule nanofilt_lite as nanofilt with:
    input:
        "results/porechop/{strain}_{stage}.fastq.gz"
    params:
        length = "8000",
        quality = get_ont_quality_filter
    output:
        "results/nanofilt/{stage}/{strain}_{stage}.fastq.gz"
    log:
        "logs/nanofilt/{strain}_{stage}.log"

"""rule coverage_ont:
    input:
        "results/reports/trimmed/nanoplot/{strain}/NanoStats.txt"
    output:
        "results/reports/trimmed/nanoplot/{strain}_coverage.txt"
    params:
        tech = "ont",
        genomesize = get_genome_size
    conda:
        "../envs/watstats.yaml"
    script:
        "../scripts/coverage_old.py"

rule coverages_ont:
    input:
        expand("results/reports/trimmed/nanoplot/{strain}_coverage.txt", strain=get_all_strain_ids())
    output:
        outfile = "results/reports/trimmed/all_reports/ont_coverages.txt"
    params:
        indir = "results/reports/trimmed/nanoplot/" #lambda wildcards, output: Path(output.outfile).parent
    shell:
        "cat {input} > {output} && "
        "rm {params.indir}*coverage.txt"
        """