#from pathlib import Path

## assembly QC
if config["checkm2_db"]["use_local"]:

    rule copy_local_checkM2_DB:
        output:
            dbfile=get_checkm2_db(),
            tar=temp(get_checkm2_tar()),
        params:
            local=config["checkm2_db"]["local_path"],
            resource_folder= get_resource_path(),
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
            dbfile=get_checkm2_db(),
        params:
            direct=get_resource_path(),
            #lambda wildcards, output: Path(output.dbfile).parent.parent,
        log:
            "logs/checkm2_DB_download.log",
        conda:
            "../envs/checkm2.yaml"
        shell:
            "checkm2 database --download --path {params.direct} > {log} 2>&1"


rule checkm2_run:
    input:
        assembly=rules.move_and_gz_assembly.output.fa_gz,
        dbfile=get_checkm2_db(),
    output:
        stats="results/{date}/qc/checkm2/{sample}/quality_report.tsv",
    params:
        outdir=lambda wildcards, output: Path(output.stats).parent,
        ext=lambda wildcards, input: Path(input.assembly).suffix,
    log:
        "logs/{date}/assembly/checkm2/{sample}.log",
    threads: 4
    conda:
        "../envs/checkm2.yaml"
    shell:
        "checkm2 predict -x {params.ext} --threads {threads} --force "
        "--database_path {input.dbfile} --input {input.assembly} "
        "--output-directory {params.outdir}/ > {log} 2>&1"
        # --remove_intermediates-x fa.gz


rule quast:
    input:
        fasta=rules.move_and_gz_assembly.output.fa,
        gff=rules.prokka.output.faa[1],
    output:
        multiext("results/{date}/qc/quast/{sample}/report.", "html", "tex", "txt", "pdf", "tsv"),
        multiext("results/{date}/qc/quast/{sample}/transposed_report.", "tex", "txt", "tsv"),
        multiext(
            "results/{date}/qc/quast/{sample}/basic_stats/",
            "cumulative_plot.pdf",
            "GC_content_plot.pdf",
        ),
        "results/{date}/qc/quast/{sample}/icarus.html",
    log:
        "logs/{date}/assembly/quast/{sample}.log",
    threads: 2
    params:
        extra="--min-contig 500 --min-identity 95.0",
    wrapper:
        "v3.3.6/bio/quast"


"""
busco -i spades/ont_66/scaffolds.fasta --auto-lineage-prok -o busco_66 -m genome
checkm lineage_wf -x fasta results/spades/ont_66/ results/checkm/ont_66/
checkm lineage_wf -x fasta results/spades/ont_66/ results/checkm/uni_156/
quast -o results/quast/66/ results/spades/66/scaffolds.fasta 

rule quast_spades:
    input:
        "results/assembly/{strain}/{strain}_scaffolds.fasta"
    output:
        repo = "results/reports/assembly/quast/{strain}_spades/report.txt"
    params:
        outdir = lambda wildcards, output: Path(output.repo).parent
    log:
        "logs/quast/spades_{strain}.log"
    conda:
        "../envs/quast.yaml"
    shell:
        "quast -o {params.outdir}/ {input} > {log} 2>&1"

rule checkm_spades:
    input:
        scaf = "results/assembly/{strain}/{strain}_scaffolds.fasta"
    output:
        repo = "results/reports/assembly/checkm/{strain}_spades/lineage.ms"
    params:
        indir = lambda wildcards, input: Path(input.scaf).parent,
        outdir = lambda wildcards, output: Path(output.repo).parent
    log:
        "logs/checkm/spades_{strain}.log"
    conda:
        "../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -x fasta {params.indir}/ {params.outdir}/ > {log} 2>&1"
rule checkm_canu:
    input:
        ass = "results/pilon/{strain}/{strain}.fasta"
    output:
        repo = "results/reports/assembly/checkm/{strain}_canu/lineage.ms"
    params:
        indir = lambda wildcards, input: Path(input.ass).parent,
        outdir = lambda wildcards, output: Path(output.repo).parent
    log:
        "logs/checkm/canu_{strain}.log"
    conda:
        "../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -x fasta {params.indir}/ {params.outdir}/ > {log} 2>&1"

rule checkm:
    input:
        ass = "results/final_assemblies/{strain}/assembly.fasta"
    output:
        repo = "results/reports/assembly/checkm/{strain}/lineage.ms"
    params:
        indir = lambda wildcards, input: Path(input.ass).parent,
        outdir = lambda wildcards, output: Path(output.repo).parent
    log:
        "logs/checkm/final_{strain}.log"
    conda:
        "../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -x fasta {params.indir} {params.outdir} 2> {log}"

rule quast_uni:
    input:
        "results/assembly/{strain}/{strain}_assembly.fasta"
    output:
        repo = "results/reports/assembly/{strain}/quast/report.tsv"
    params:
        outdir = lambda wildcards, output: Path(output.repo).parent,
        busco = "-b"#"--conserved-genes-finding"
    log:
        "logs/quast/unicycler_{strain}.log" # log file automated created - need to mv it
    conda:
        "../envs/quast.yaml"
    shell:
        "quast {params.busco} -o {params.outdir}/ {input} 2> {log}"

rule checkm_uni:
    input:
        ass = "results/assembly/{strain}/{strain}_assembly.fasta"
    output:
        repo = "results/reports/assembly/{strain}/checkm/lineage.ms"
    params:
        indir = lambda wildcards, input: Path(input.ass).parent,
        outdir = lambda wildcards, output: Path(output.repo).parent
    log:
        "logs/checkm/uni_{strain}.log"
    conda:
        "../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -x fasta {params.indir}/ {params.outdir}/ 2> {log}"
"""
'''
rule busco:
    input:
        "results/assembly/{strain}/assembly.fasta"
    output:
        dir("results/reports/assembly/{strain}/busco/")
    params:
        outfolder_name = "{strain}_busco"
    #logfile is produced in log subfolder -> need to mv it
    shell:
        "busco -i {input} -o {params.outdir} -m genome --auto-lineage-prok "
        "&& mv {params.outfolder_name}/* {output}/*"
'''