from pathlib import Path


rule prokka:
    input:
        rules.unicycler.output.assembly,
    output:
        faa="results/{date}/analysis/prokka/{sample}/{sample}.faa",
    params:
        outdir=lambda wildcards, output: Path(output.faa).parent,
    log:
        "logs/{date}/analysis/prokka/{sample}.log",
    conda:
        "../envs/prokka.yaml"
    shell:
        "prokka --outdir {params.outdir}/ --force "
        "--prefix {wildcards.sample} {input} > {log} 2>&1"  #--quiet


# for metagenomic data:
# --metagenome


rule clone_plm_arg:
    output:
        main=get_plm_arg_main(),
        dummy=touch("clone_plm_arg.done"),
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
            model=get_plm_arg_model_file(),
        params:
            local=config["plm_arg"]["model"]["local_path"],
            folder=get_plm_arg_model_path(),
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
            model=get_plm_arg_model_file(),
        params:
            download=config["plm_arg"]["model"]["url"],
            folder=get_plm_arg_model_path(),
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
            folder=get_plm_arg_model_path(),
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
            reg=get_plm_arg_regression_file(),
        params:
            download=config["plm_arg"]["regression"]["url"],
            folder=get_plm_arg_model_path(),
        log:
            "logs/plm_arg_regression_download.log",
        conda:
            "../envs/unix.yaml"
        shell:
            "(cd {params.folder} && "
            "wget {params.download}) > {log} 2>&1"


rule run_plm_arg:
    input:
        model=get_plm_arg_model_file(),
        reg=get_plm_arg_regression_file(),
        main=get_plm_arg_main(),
        fasta=rules.prokka.output.faa,
    output:
        "results/{date}/analysis/plm_arg/{sample}.tsv",
    params:
        folder=lambda wildcards, input: Path(input.main).parent,
        main= lambda wildcards, input: Path(input.main).name,
        root=get_root(),
    log:
        "logs/{date}/analysis/plm_arg/{sample}.log",
    threads: 10
    conda:
        "../envs/plm_arg.yaml"
    shell:
        "(cd {params.folder} && "
        "python {params.main} predict -i {params.root}/{input.fasta} "
        "-o {params.root}/{output} > {log} 2>&1) "


"""
rule rename_assembly:
    input:
        "results/final_assemblies/{strain}/assembly.fasta"
    output:
        "results/assembly/{strain}_contig.fasta"
    shell:
        "cp {input} {output}"

rule prokka:
    input:
        "results/assembly/{strain}/{strain}_assembly.fasta"
        #"results/assembly/{strain}_contig.fasta"
    output:
        gff = "results/analysis/prokka/{strain}/{strain}.gff"
    params:
        outdir = lambda wildcards, output: Path(output.gff).parent
    log:
        "logs/prokka/{strain}.log"
    conda:
        "../envs/prokka.yaml"
    shell:
        "prokka --outdir {params.outdir}/ --force --prefix {wildcards.strain} {input} 2> {log}"

rule abricate:
    input:
        "results/assembly/{strain}/{strain}_assembly.fasta"
        #"results/assembly/{strain}_contig.fasta"
    output:
        all = "results/analysis/abricate/{strain}/{strain}_all.csv",
        info = "results/analysis/abricate/{strain}/{strain}_info.txt"
        #ncbi = "results/analysis/abricate/{strain}/{strain}_ncbi.tab",
        #summary = "results/analysis/abricate/{strain}/{strain}_summary.txt"
    params:
        outdir = lambda wildcards, output: Path(output.all).parent
    log:
        "logs/abricate/{strain}.log"
    conda:
        "../envs/abricate.yaml"
    script:
        "../scripts/abricate_all.py"
        
## visualize the annotated genome
rule gview:
    input:
        gff = "results/analysis/prokka/{strain}/{strain}.gff",
        gbk = "results/analysis/prokka/{strain}/{strain}.gbk"
    output:
        "results/reports/assembly/{strain}_gview.png"
    params:
        extra = "-l circular",
        jar_file = "resources/gview/gview.jar",
        ## style sheet
        gss = "resources/gview/example_styles/gssExample.gss"
        ## we need a style sheet to make this look cool! (see basic and gss examples)
        #style_sheet = 
    log:
        "logs/gview/{strain}.log"
    conda:
        "../envs/gview.yaml" # env with java
    shell:
        "java -jar {params.jar_file} -i {input.gbk} -o {output} "
        "-g {input.gff} -s {params.gss} {params.extra}"
  """
