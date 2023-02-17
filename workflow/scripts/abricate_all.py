'''TODO:
include log file
abricate db from config file
'''

import os
import pandas as pd

def make_all_file(outdir, all_path, info_out_path, strain, dbs, assembly, logfile):
    first_file = True
    info_ls = []

    ## runs abricate for each database in the DB list & saves all information 
    for db in dbs:
        abricate_out = "{}/{}_{}.tsv".format(outdir, strain, db)
        abricate(db, assembly, abricate_out, logfile)

        tsv_df = pd.read_csv(abricate_out, sep="\t")
        info_ls.append((db, len(tsv_df)))
        
        if not tsv_df.empty:
            if first_file:
                write_info_file(info_ls, info_out_path, assembly, first_file)

                all_df = tsv_df.copy()
                first_file = False
            else:
                all_df = pd.concat([all_df,tsv_df], ignore_index=True)
        ## after saving the information the abricate output file is deleted 
        remove_abricate(abricate_out)

    ## write output files
    write_info_file(info_ls, info_out_path, assembly, first_file)
    all_df.to_csv(all_path)

def abricate(db, assembly, outpath, logfile):
    os.system(f"abricate --db {db} {assembly} > {outpath} 2>> {logfile}") #funktionierts mit 2 >?

def write_info_file(info_ls, info_out_path, assembly, first):
    if first:
        info_file = open(info_out_path, "w")
        info_file.write("file: {}\n\n".format(assembly))
        info_file.close()
    else:
        info_df = pd.DataFrame.from_records(info_ls, columns=["database", "#results"])
        info_df.to_csv(info_out_path, mode="a", sep="\t", index=False)

def remove_abricate(file):
    os.system("rm {}".format(file))

'''
def write_summary(outdir, strain):
    abricate_outfiles = glob.glob(outdir + "/*.tab")
    gene_dict = {}
    for outfile in abricate_outfiles:
        tabfile = open(outfile, "r")
        for line in tabfile:
            if line.startswith("#"):
                continue

            gene_start = line.split()[2]
            if gene_start in gene_dict:
                gene_dict[gene_start].append(line)
            else:
                gene_dict[gene_start] = [line]
        tabfile.close()

    summary =f"{outdir}/{strain}_summary.txt"
    sumfile = open(summary, "w")
    for gene in gene_dict:
        found_in_dbs = len(gene_dict[gene])
        sumfile.write(f"genes starting at position {gene} found in {found_in_dbs}" + "\n")
        for line in gene_dict[gene]:
            sumfile.write(line)

    sumfile.close()'''

dbs = ["ncbi", "plasmidfinder", "argannot", "megares", "ecoli_vf", "resfinder", "card", "ecoh", "vfdb"]


assembly = str(snakemake.input) #"results/final_assemblies/115_copy/assembly.fasta" #str(snakemake.input)
all_path = str(snakemake.output.all) #"{}/{}_all.csv".format(outdir, strain)
info_out_path = str(snakemake.output.info) #"{}/{}_info.txt".format(outdir, strain)
outdir = str(snakemake.params.outdir) #"results/analysis/abricate/115/" #str(snakemake.params.outdir)
strain = str(snakemake.wildcards.strain)
logfile = str(snakemake.log)


make_all_file(outdir, all_path, info_out_path, strain, dbs, assembly, logfile)
