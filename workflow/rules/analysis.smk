from pathlib import Path


rule prokka:
    input:
        local(get_assembly),
    output:
        faa=local(
            multiext("results/{date}/analysis/prokka/{sample}/{sample}.", "faa", "gff")
        ),
    params:
        outdir=lambda wildcards, output: Path(output.faa[0]).parent,
    log:
        "logs/{date}/analysis/prokka/{sample}.log",
    threads: 64
    conda:
        "../envs/prokka.yaml"
    shell:
        "prokka --outdir {params.outdir}/ --force "
        "--prefix {wildcards.sample} --cpus {threads} "
        "{input} > {log} 2"
        #--quiet


rule load_genomad_DB:
    output:
        file=local(get_genomad_DB_file()),
    params:
        folder=lambda wildcards, output: Path(output.file).parent.parent,
    log:
        "logs/load_genomad_DB.log",
    conda:
        "../envs/genomad.yaml"
    shell:
        "genomad download-database {params.folder}/ > {log} 2>&1"


rule genomad_run:
    input:
        db=local(rules.load_genomad_DB.output.file),
        asmbl=rules.assembly_gz.output.fa_gz,
    output:
        plasmid_tsv="results/{date}/analysis/genomad/{sample}/{sample}_summary/{sample}_plasmid_summary.tsv",
        virus_tsv="results/{date}/analysis/genomad/{sample}/{sample}_summary/{sample}_virus_summary.tsv",
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
        "{input.asmbl} {params.outdir}/ "
        "{params.db_folder}/ > {log} 2>&1"


# for metagenomic data:
# --metagenome
if config["card"]["data"]["use_local"]:

    rule copy_local_CARD_data:
        output:
            json=get_card_db_file(),
        params:
            local=config["card"]["data"]["local_path"],
            folder=local(lambda wildcards, output: Path(output.json).parent),
            tar_name=get_card_tar_file(),
        log:
            "logs/CARD_data_local_copy.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(mkdir -p {params.folder}/ && "
            "cd {params.folder}/ && "
            "cp {params.local} . && "
            "tar -xvf {params.tar_name}) > {log} 2>&1"


"""
else:

    rule download_CARD_data:
        output:
            model=get_card_db_folder(),
        params:
            download=config["card"]["data"]["url"],
            folder=get_card_db_folder(),
        log:
            "logs/CARD_data_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"
"""

# rgi load -i resources/card_db/card.json --local
# rgi main -i results/test_115/analysis/prokka/ab_115/ab_115.faa -o results/test_115/analysis/rgi_115.out -t protein -a DIAMOND -d wgs --clean --local


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
        "rgi load -i {input} --local > {log} 2>&1"  # 


rule CARD_run:
    input:
        faa=local(rules.prokka.output.faa[0]),
        db=local(rules.CARD_load_DB.output),
    output:
        txt="results/{date}/analysis/card/{sample}.txt",
    params:
        path_wo_ext=lambda wildcards, output: local(Path(output.txt).with_suffix("")),
    log:
        "logs/{date}/analysis/card/{sample}.log",
    threads: 64
    conda:
        "../envs/card.yaml"
    shell:
        "rgi main -i {input.faa} -o {params.path_wo_ext} -t protein "
        "-a DIAMOND -d wgs --local --clean > {log} 2>&1"
        # 


rule clone_plm_arg:
    output:
        main=local(get_plm_arg_main()),
        dummy=touch("logs/clone_plm_arg.done"),
    log:
        "logs/plm_arg_clone.log",
    params:
        git_url=config["plm_arg"]["git"],
        folder=lambda wildcards, output: Path(output.main).parent,
    conda:
        "../envs/unix.yaml"
    shell:
        "(rm -r -f {params.folder} && "
        "git clone --recursive {params.git_url} {params.folder}) > {log} 2>&1"


if config["plm_arg"]["model"]["use_local"]:

    rule copy_local_plm_arg_model:
        input:
            rules.clone_plm_arg.output.dummy,
        output:
            model=local(get_plm_arg_model_file()),
        params:
            local=config["plm_arg"]["model"]["local_path"],
            folder=lambda wildcards, output: Path(output.model).parent,
        log:
            "logs/plm_arg_model_local_copy.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "cp {params.local} {params.folder}/ > {log} 2>&1"

else:

    rule download_plm_arg_model:
        input:
            rules.clone_plm_arg.output.dummy,
        output:
            model=local(get_plm_arg_model_file()),
        params:
            download=config["plm_arg"]["model"]["url"],
            folder=lambda wildcards, output: Path(output.model).parent,
        log:
            "logs/plm_arg_model_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


if config["plm_arg"]["regression"]["use_local"]:

    rule copy_local_plm_arg_regression:
        input:
            rules.clone_plm_arg.output.dummy,
        output:
            reg=get_plm_arg_regression_file(),
        params:
            local=config["plm_arg"]["regression"]["local_path"],
            folder=lambda wildcards, output: Path(output.reg).parent,
        log:
            "logs/plm_arg_regression_local_copy.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "cp {params.local} {params.folder}/ > {log} 2>&1"

else:

    rule download_plm_arg_regression:
        input:
            rules.clone_plm_arg.output.dummy,
        output:
            reg=local(get_plm_arg_regression_file()),
        params:
            download=config["plm_arg"]["regression"]["url"],
            folder=lambda wildcards, output: Path(output.reg).parent,
        log:
            "logs/plm_arg_regression_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


rule run_plm_arg:
    input:
        model=local(get_plm_arg_model_file()),
        reg=local(get_plm_arg_regression_file()),
        main=local(get_plm_arg_main()),
        fasta=rules.prokka.output.faa[0],
    output:
        tsv=local("results/{date}/analysis/plm_arg/{sample}/{sample}.tsv"),
    params:
        folder=lambda wildcards, input: Path(input.main).parent,
        main=lambda wildcards, input: Path(input.main).name,
        root=local(get_root()),
    log:
        local("logs/{date}/analysis/plm_arg/{sample}.log"),
    threads: 30
    conda:
        "../envs/plm_arg.yaml"
    shell:
        "(cd {params.folder}/ && "
        "python {params.main} predict -i {params.root}/{input.fasta} "
        "-o {params.root}/{output.tsv}) > {log} 2>&1"


rule extract_plm_arg_results:
    input:
        tsv=local(rules.run_plm_arg.output.tsv),
    output:
        arg="results/{date}/analysis/plm_arg/{sample}/{sample}_arg.csv",
        non_arg="results/{date}/analysis/plm_arg/{sample}/{sample}_non_arg.csv",
    log:
        "logs/{date}/analysis/plm_arg/{sample}_extract.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/extract_plm_arg.py"


rule plot_plm_arg_results:
    input:
        csv=expand(
            "results/{{date}}/analysis/plm_arg/{sample}/{sample}_arg.csv",
            sample=get_samples(),
        ),
    output:
        html="results/{date}/report/arg_plot.html",
    log:
        "logs/{date}/report/arg_plot.log",
    threads: 2
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/plot_resistance.py"
