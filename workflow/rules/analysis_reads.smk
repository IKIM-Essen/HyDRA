## convert fastq to fasta
rule fastq_to_fasta:
    input:
        fastq=rules.chopper.output.trim_filt,
    output:
        fasta = "results/{date}/trimmed/chopper/{sample}.fasta",
    log:
        "logs/{date}/ont/fastq_to_fasta/{sample}.log",
    threads: 32
    conda:
        "../envs/unix.yaml"
    shell:
        "zcat {input.fastq} | sed -n '1~4s/^@/>/p;2~4p' > {output.fasta} 2> {log}"


rule load_genomad_DB:
    output:
        file=get_genomad_DB_file(),
    params:
        folder=lambda wildcards, output: Path(output.file).parent.parent,
    log:
        "logs/load_genomad_DB.log",
    conda:
        "../envs/genomad.yaml"
    shell:
        "genomad download-database {params.folder}/ > {log} 2>&1"


rule genomad_reads_run:
    input:
        db=rules.load_genomad_DB.output.file,
        reads=rules.fastq_to_fasta.output.fasta,
    output:
        plasmid_tsv="results/{date}/analysis/plasmids/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv",
        virus_tsv="results/{date}/analysis/plasmids/{sample}/{sample}_summary/{sample}_virus_summary.tsv",
    params:
        db_folder=lambda wildcards, input: Path(input.db).parent,
        outdir=lambda wildcards, output: Path(output.plasmid_tsv).parent.parent,
    log:
        "logs/{date}/analysis/genomad/{sample}.log",
    threads: 64
    conda:
        "../envs/genomad.yaml"
    shell:
        "genomad end-to-end --cleanup -t {threads} "
        "{input.reads} {params.outdir}/ "
        "{params.db_folder}/ > {log} 2>&1"


rule download_CARD_data:
        output:
            json = get_card_db_file(),
        params:
            download = config["card"]["data"]["url"],
            folder = lambda wildcards, output: Path(output.json).parent,
            filename = lambda wildcards, output: Path(output.json).name,
        log:
            "logs/CARD_data_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cd {params.folder} && "
            "wget {params.download} && "
            "tar -xvf data ./{params.filename}) > {log} 2>&1"


rule CARD_load_DB:
    input:
        get_card_db_file(),
    output:
        touch("logs/CARD_load_DB.done"),
    log:
        "logs/CARD_load_DB.log",
    conda:
        "../envs/card.yaml"
    shell:
        "rgi clean --local && "
        "rgi load --card_json {input} --local > {log} 2>&1"


rule CARD_annotation:
    input:
        json = get_card_db_file(),
    output:
        touch("logs/CARD_annotation.done"),
    params:
        folder = lambda wildcards, input: Path(input.json).parent,
        file = lambda wildcards, input: Path(input.json).name,
        ann = get_card_annotation_file(),
    log:
        "logs/CARD_annotation.log",
    conda:
        "../envs/card.yaml"
    shell:
        "(cd {params.folder}/ && "
        "rgi card_annotation -i {params.file}) && "
        "rgi load -i {input.json} --card_annotation {params.ann} --local > {log} 2>&1"


rule CARD_reads_run:
    input:
        fastq=rules.chopper.output.trim_filt,
        db=rules.CARD_annotation.output
    output:
        txt="results/{date}/analysis/resistance/{sample}/{sample}.allele_mapping_data.txt",
    params:
        folder=lambda wildcards, output: Path(output.txt).parent,
    log:
        "logs/{date}/analysis/card/{sample}.log",
    threads: 30
    conda:
        "../envs/card.yaml"
    shell:
        "rgi bwt -1 {input.fastq} "
        "-o {params.folder}/{wildcards.sample} --aligner bwa "
        "--clean -n {threads} --local > {log} 2>&1"
