import pandas as pd
import glob

configfile: "config/config.yaml"

samples = pd.read_csv("config/pep/documents.csv").set_index("sample_name", drop=False)
samples.index.names = ["sample_name"]
ill_data_path = config["ill_data_path"]
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

def get_nanopore_barcode(wildcards):
    return pep.sample_table.loc[wildcards.strain][["nanopore_barcode"]]

def get_genome_size(wildcards):
    return pep.sample_table.loc[wildcards.strain]["species_genome_size"]

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
        read_path = ill_data_path + "{strain}_{lane}_{read}.fastq.gz"
        return read_path
    elif wildcards.stage == stages[1]:
        return  "results/preprocess_ill/{strain}/{strain}_{read}_trimmed.fastq.gz"

def get_ill_rawR1(wildcards):
    read_path = ill_data_path + "{strain}_{lane}_R1.fastq.gz"
    return read_path

def get_ill_rawR2(wildcards):
    read_path = ill_data_path + "{strain}_{lane}_R2.fastq.gz"
    return read_path

def get_summaryfile_by_run(wildcards):
    if wildcards.run == ont_runs[0]:
        return "/projects/seqlab/GridIONs/20210713_JW001_KID1-3_Ec_ZKB_NBD104/20210713_JW001_KID1-3_Ec_ZKB_NBD104/20210713_1530_X2_FAP91220_d47f932e/sequencing_summary_FAP91220_e888cbe6.txt"
    elif wildcards.run == ont_runs[1]:
        return "/projects/seqlab/GridIONs/20211115_SM0006_JW0002_CefRes_strains_20211115/SM0006_JW0002_CefRes_strains_20211115/20211115_1648_X1_FAP91415_d0c3f758/sequencing_summary_FAP91415_f2d1671c.txt"

"""
if wildcards.lane == lanes[0]:
        return pep.sample_table.loc[wildcards.strain][["illumina_read2l1"]]
    elif wildcards.lane == lanes[1]:
        return pep.sample_table.loc[wildcards.strain][["illumina_read2l2"]]
"""   