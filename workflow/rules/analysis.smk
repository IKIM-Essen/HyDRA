from pathlib import Path

rule rename_assembly:
    input:
        "results/final_assemblies/{strain}/assembly.fasta"
    output:
        "results/assembly/{strain}_contig.fasta"
    shell:
        "cp {input} {output}"

rule prokka:
    input:
        "results/assembly/{strain}_contig.fasta"
    output:
        gff = "results/analysis/prokka/{strain}/{strain}.gff"
    params:
        outdir = lambda wildcards, output: Path(output.gff).parent
    log:
        "logs/prokka/{strain}.log"
    conda:
        "../envs/prokka.yaml"
    shell:
        "prokka --outdir {params.outdir}/ --force --prefix {wildcards.strain} {input} > {log} 2>&1"

rule abricate:
    input:
        "results/assembly/{strain}_contig.fasta"
    output:
        ncbi = "results/analysis/abricate/{strain}/{strain}_ncbi.tab",
        summary = "results/analysis/abricate/{strain}/{strain}_summary.txt"
    params:
        outdir = lambda wildcards, output: Path(output.ncbi).parent
    log:
        "logs/abricate/{strain}.log"
    conda:
        "../envs/abricate.yaml"
    script:
        "../scripts/abricate_summary.py"