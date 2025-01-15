import os


configfile: "config/config.yaml"


## helper functions
def get_root():
    return os.getcwd()


def get_assembly_type():
    return config["assembly_type"]


def get_has_short_reads():
    if get_assembly_type() in ["hybrid", "short"]:
        return True
    else:
        return False


def get_has_long_reads():
    if get_assembly_type() in ["hybrid", "long"]:
        return True
    else:
        return False


def get_data_path_ill():
    return config["data_handling"]["data"]["illumina"]


def get_data_path_ont():
    return config["data_handling"]["data"]["ont"]


def get_resource_path():
    return config["data_handling"]["resources"]


def get_run_date():
    return config["run_date"]


def get_samples():
    return list(pep.sample_table["sample_name"].values)


def get_ill_fastqs(wildcards):
    return (
        pep.sample_table.loc[wildcards.sample]["shortR1"],
        pep.sample_table.loc[wildcards.sample]["shortR2"],
    )


def get_ont_fastq(wildcards):
    return pep.sample_table.loc[wildcards.sample]["long"]


def get_adapters():
    return config["adapter_seqs"]


def get_multiqc_config():
    return "config/multiqc_config.yaml"


def get_multiqc_input(wildcards):
    if get_has_short_reads():
        short_in = expand(
            [
                "results/{{date}}/trimmed/fastp/{sample}.fastp.json",
                "results/{{date}}/qc/fastqc/{sample}/ill_{sample}_{read}_fastqc.zip",
            ],
            sample=get_samples(),
            read=["1", "2"],
        )

    if get_has_long_reads():
        long_in = expand(
            [
                "results/{{date}}/qc/nanoplot/{sample}/{sample}_NanoStats.txt",
            ],
            sample=get_samples(),
        )
    if get_assembly_type() == "hybrid":
        short_in.extend(long_in)
        return short_in
    elif get_assembly_type() == "short":
        return short_in
    elif get_assembly_type() == "long":
        return long_in


def get_assembly(wildcards):
    return "results/{date}/assembly/{sample}/assembly.fasta"


def get_checkm2_db():
    path = "{0}{1}".format(get_resource_path(), config["checkm2_db"]["db_file"])
    return path


def get_checkm2_tar():
    name = (config["checkm2_db"]["local_path"]).split("/")[-1]
    path = "{0}{1}".format(get_resource_path(), name)
    return path


def get_genomad_DB_file():
    path = "{}genomad_db/names.dmp".format(get_resource_path())
    return path


def get_card_db_file():
    name = config["card"]["data"]["dbfile"]
    path = "{}CARD_db/{}".format(get_resource_path(), name)
    return path


def get_card_tar_file():
    if config["card"]["data"]["use_local"]:
        name = Path(config["card"]["data"]["local_path"]).name
    else:
        name = "card-data.tar.bz2"
    # path = "{}CARD_db/{}".format(get_resource_path(), name)
    return name


def get_plm_arg_main():
    script = config["plm_arg"]["main"]
    script_path = "{}{}".format(get_resource_path(), script)
    return script_path


def get_plm_arg_model_path():
    path = "{}PLM-ARG/models/".format(get_resource_path())
    return path


def get_plm_arg_model_file():
    name = (config["plm_arg"]["model"]["url"]).split("/")[-1]
    path = "{0}{1}".format(get_plm_arg_model_path(), name)
    return path


def get_plm_arg_regression_file():
    name = (config["plm_arg"]["regression"]["url"]).split("/")[-1]
    path = "{0}{1}".format(get_plm_arg_model_path(), name)
    return path
