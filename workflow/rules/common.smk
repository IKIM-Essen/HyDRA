import pandas as pd
import glob

configfile: "config/config.yaml"

samples = pd.read_csv("config/pep/documents.csv").set_index("sample_name", drop=False)
samples.index.names = ["sample_name"]
stages = ["before_trim", "trimmed"]
lanes = ["L001","L002"]
read_ids = ["R1","R2"]
ont_runs = ["run_1", "run_2"]
filter_less = ["109", "102", "265", "188", "167", "156", "139", "129"]
filter_more = ["115", "120"]

## helper functions

def get_all_strain_ids():
    return pep.sample_table["sample_name"].to_list()

def get_nanopore_reads(wildcards):
    return pep.sample_table.loc[wildcards.strain][["nanopore_reads"]]

'''def get_nanopore_barcode(wildcards):
    return pep.sample_table.loc[wildcards.strain][["nanopore_barcode"]]
'''
def get_genome_size(strain):
    return pep.sample_table.loc[strain]["species_genome_size"]

def get_genome_size_dict():
    genome_size_dict = {}
    for sample in get_all_strain_ids():
        genome_size_dict[sample] = get_genome_size(sample)
    return genome_size_dict

def get_all_stages():
    return stages

def get_all_lanes():
    return lanes

def get_all_ont_runs():
    return ont_runs

def get_all_read_ids():
    return read_ids

def get_ont_quality_filter(wildcards):
    if wildcards.strain in filter_less:
        quality = "14"
        return quality
    elif wildcards.strain in filter_more:
        quality = "16"
        return quality
    else:
        quality = "15"
        return quality

def get_ont_reads_by_stage(wildcards):
    if wildcards.stage == stages[0]:
        return "data/ont/{strain}.fastq.gz"
    elif wildcards.stage == stages[1]:
        return "results/nanofilt/{stage}/{strain}_{stage}.fastq.gz"

def get_ill_reads_by_stage(wildcards):
    if wildcards.stage == stages[0]:
        read_path = "data/illumina/{strain}_{lane}_{read}.fastq.gz"
        return read_path
    elif wildcards.stage == stages[1]:
        return  "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq.gz"

def get_ill_rawR1(wildcards):
    read_path = "data/illumina/{strain}_{lane}_R1.fastq.gz"
    return read_path

def get_ill_rawR2(wildcards):
    read_path = "data/illumina/{strain}_{lane}_R2.fastq.gz"
    return read_path

def get_summaryfile_by_run(wildcards):
    if wildcards.run == ont_runs[0]:
        return config["summary_ont_run1"]
    elif wildcards.run == ont_runs[1]:
        return config["summary_ont_run2"]
"""
if wildcards.lane == lanes[0]:
        return pep.sample_table.loc[wildcards.strain][["illumina_read2l1"]]
    elif wildcards.lane == lanes[1]:
        return pep.sample_table.loc[wildcards.strain][["illumina_read2l2"]]
"""   