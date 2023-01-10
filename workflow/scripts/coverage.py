import pandas as pd
import subprocess

''' reads total number of reads & the average read length from the multiqc fatsqc file
    renames the columns & saves it to a DataFrame '''
def make_ill_df(multiqc_fastqc):

    ## relevant columns for Illumina reads are loaded to ill_df & renamed
    ill_cols = ["Sample", "avg_sequence_length", "Total Sequences"]
    ill_col_rename = {ill_cols[1]: "Ill_avg_read_len", ill_cols[2]:"Ill_total_reads"}
    
    ill_df = pd.read_table(multiqc_fastqc, index_col="Sample", usecols=ill_cols).rename(columns=ill_col_rename)

    ## since there are 2 fastqc files per sample (R1 & R2), need to remove the R2 rows
    new_ind = {}
    rm_lst = []
    ## if stage ="trimmed" just R1 is used
    if stage == "trimmed":
        for ind in ill_df.index.tolist():
            if ind.find("R1")>=0:
                new_ind[ind] = ind[:ind.rfind("_R1")]
            else:
                rm_lst.append(ind)

        ill_df = ill_df.drop(index=rm_lst).rename(index=new_ind)

    ## if stage = before_trim there are 2 lanes per R1, need to sum up the total reads of L001 & L002
    elif stage == "before_trim":
        new_data = []
        for ind in ill_df.index.tolist():
    
            if ind.find("R2")>=0:
                rm_lst.append(ind)
            else:
                if ind.find("L001")>=0:
                    total_reads = ill_df["Ill_total_reads"].loc[ind]
                    new_ind[ind] = ind[:ind.rfind("_L001")]
                else:
                    total_reads += ill_df["Ill_total_reads"].loc[ind]
                    new_data.append(total_reads)
                    rm_lst.append(ind)

        ill_df = ill_df.drop(index=rm_lst).rename(index=new_ind)
        ill_df["Ill_total_reads"] = new_data

    return(ill_df)

''' reads total number of reads & the average read length from the multiqc nanostats file
    renames the columns & saves it to a DataFrame '''
def make_ont_df(multiqc_nanostats):

    ## relevant columns for ONT reads are loaded to ont_df & renamed
    ont_cols = ["Sample", "Mean read length_fastq", "Number of reads_fastq"]
    ont_col_rename = {ont_cols[1]:"ONT_avg_read_len", ont_cols[2]:"ONT_total_reads"}

    ont_df = pd.read_table(multiqc_nanostats, index_col="Sample", usecols=ont_cols).rename(columns=ont_col_rename)

    ## since sample names are numbers, need to be converted to str
    if stage == "trimmed":
        ont_df.index = ont_df.index.astype(str)

    elif stage == "before_trim":
        new_ind = {}
        for ind in ont_df.index.tolist():
            new_ind[ind] = ind[:ind.rfind("_before")]
        ont_df = ont_df.rename(index=new_ind)

    return(ont_df)

def append_to_dict_lst(cov_dict, col, value):
    ls = cov_dict[col]
    ls.append(value)
    cov_dict[col] = ls
    return(cov_dict)

''' calculates coverage with watstats, reads coverage & expected percentage genome covered from STDOUT
    & saves them with the samplename and the genome size to a dataframe ''' 
def coverage_to_df(df, watstats_path, techs, genome_sizes):
    cov_dict = {}
    first_sample = True

    for sample in df.index.tolist():
        genome_size = str(genome_sizes[sample])

        if first_sample:
            cov_dict["Sample"] = [sample]
            cov_dict["Genome_size"] = [genome_size]
        else:
            cov_dict = append_to_dict_lst(cov_dict, "Sample", sample)
            cov_dict = append_to_dict_lst(cov_dict, "Genome_size", genome_size)

        for tech in techs:
            ## get the values for watstats calculation
            num_reads = str(df.at[sample,f"{tech}_total_reads"])
            avg_len = str(df.at[sample,f"{tech}_avg_read_len"])

            #watstats cmd: perl {watstats_path} -g {genome_size} -n {num_reads} -l {avg_len}
            cmd = ["perl", watstats_path, "-g", genome_size, "-n", num_reads, "-l", avg_len]
            out_txt = subprocess.check_output(cmd).decode()

            ## extract values from watstats output & save in dict (key = samplename)
            out_txt = out_txt.strip().split('\n')
            for elem in out_txt:
                if elem.find("coverage")>=0:
                    elem_split = elem.split()
                    value = elem_split[-1]
                    if "Percentage" in elem_split:
                        name = tech + "_percentage_genome_covered"
                    elif "Redundancy" in elem_split:
                        name = tech + "_coverage"
                    
                    if first_sample:
                        cov_dict[name] = [value]
                    else:
                        cov_dict = append_to_dict_lst(cov_dict, name, value)
        first_sample = False

    ## save the dict to a dataframe and return it
    cov_df = pd.DataFrame.from_dict(cov_dict)
    cov_df.set_index("Sample", inplace=True)
    return(cov_df)

# reorder the columns in the dataframe before saving it to csv
def reorder_output(all_df):
    col_order = ["_total_reads", "_avg_read_len", "_coverage", "_percentage_genome_covered"]
    all_cols = ["Genome_size"]
    for tech in techs:
        for col in col_order:
            all_cols.append(tech + col)

    return(all_df[all_cols])

watstats_path = "resources/watstats"
techs = ["ONT", "Ill"]

genome_sizes = snakemake.params.genomesizes
stage = snakemake.wildcards.stage
multiqc_nanostats = snakemake.input.ont 
multiqc_fastqc = snakemake.input.ill 
outfile = str(snakemake.output) 

ont_df = make_ont_df(multiqc_nanostats)
ill_df = make_ill_df(multiqc_fastqc)

## join ONT and Ill dfs
df = ont_df.join(ill_df)
## all values without decimal point
df = df.round(0).astype(int)

cov_df = coverage_to_df(df, watstats_path, techs, genome_sizes)
all_df = df.join(cov_df)

out_df = reorder_output(all_df)
## final output: genome size, ONT & then Ill: total reads, avg read length, coverage, & percentage covered genome
pd.DataFrame.to_csv(out_df, outfile)
