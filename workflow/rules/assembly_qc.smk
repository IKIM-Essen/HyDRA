# from pathlib import Path

## assembly QC
if config["checkm2_db"]["use_local"]:

    rule copy_local_checkM2_DB:
        output:
            dbfile=local(get_checkm2_db()),
            tar=local(get_checkm2_tar()),
        params:
            local=config["checkm2_db"]["local_path"],
            resource_folder=lambda wildcards, output: Path(output.dbfile).parent.parent,
        log:
            "logs/checkm2_DB_local_copy.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cp {params.local} {params.resource_folder}/ && "
            "tar fzxv {output.tar} -C {params.resource_folder}) > {log} 2>&1"

else:

    rule checkm2_DB_download:
        output:
            dbfile=local(get_checkm2_db()),
        params:
            direct=lambda wildcards, output: Path(output.dbfile).parent.parent,
        log:
            "logs/checkm2_DB_download.log",
        conda:
            "../envs/checkm2.yaml"
        shell:
            "checkm2 database --download --path {params.direct} > {log} 2>&1"


rule checkm2_run:
    input:
        assembly=rules.assembly_gz.output.fa_gz,
        dbfile=local(get_checkm2_db()),
    output:
        stats="results/{date}/qc/checkm2/{sample}/quality_report.tsv",
    params:
        outdir=lambda wildcards, output: Path(output.stats).parent,
        ext=lambda wildcards, input: Path(input.assembly).suffix,
    log:
        "logs/{date}/assembly/checkm2/{sample}.log",
    threads: 30
    conda:
        "../envs/checkm2.yaml"
    shell:
        "checkm2 predict -x {params.ext} --threads {threads} --force "
        "--database_path {input.dbfile} --input {input.assembly} "
        "--output-directory {params.outdir}/ > {log} 2>&1"
        # --remove_intermediates-x fa.gz


rule quast:
    input:
        fasta=get_assembly,
        gff=local(rules.prokka.output.faa[1]),
    output:
        multiext(
            "results/{date}/qc/quast/{sample}/report.",
            "html",
            "tex",
            "txt",
            "pdf",
            "tsv",
        ),
        multiext(
            "results/{date}/qc/quast/{sample}/transposed_report.", "tex", "txt", "tsv"
        ),
        multiext(
            "results/{date}/qc/quast/{sample}/basic_stats/",
            "cumulative_plot.pdf",
            "GC_content_plot.pdf",
        ),
        "results/{date}/qc/quast/{sample}/icarus.html",
    log:
        "logs/{date}/assembly/quast/{sample}.log",
    threads: 20
    params:
        extra="--min-contig 300 --min-identity 95.0",
    wrapper:
        "v3.10.2/bio/quast"
