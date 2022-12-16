from pathlib import Path
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

"""

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
        "checkm lineage_wf -x fasta {params.indir} {params.outdir} > {log} 2>&1"

rule quast_uni:
    input:
        "results/unicycler/{strain}/assembly.fasta"
    output:
        repo = "results/reports/assembly/quast/{strain}_uni/report.txt"
    params:
        outdir = lambda wildcards, output: Path(output.repo).parent
    log:
        "logs/quast/unicycler_{strain}.log"
    conda:
        "../envs/quast.yaml"
    shell:
        "quast -o {params.outdir}/ {input} > {log} 2>&1"

rule checkm_uni:
    input:
        ass = "results/unicycler/{strain}/assembly.fasta"
    output:
        repo = "results/reports/assembly/checkm/{strain}_uni/lineage.ms"
    params:
        indir = lambda wildcards, input: Path(input.ass).parent,
        outdir = lambda wildcards, output: Path(output.repo).parent
    log:
        "logs/checkm/uni_{strain}.log"
    conda:
        "../envs/checkm.yaml"
    shell:
        "checkm lineage_wf -x fasta {params.indir}/ {params.outdir}/ > {log} 2>&1"
