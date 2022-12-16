import glob
import os
import pandas as pd

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

    sumfile.close()

def combined_resultFiles(inpath):
    abricate_files = glob.glob(inpath)
    
    first = True
    for abricate_file in abricate_files:
        if first:
            all_df = pd.read_table(abricate_file, header=0)
            first = False
        else:
            db_df = pd.read_table(abricate_file, header=0)
            all_df = pd.concat([all_df, db_df], ignore_index=True)
    return(all_df)

def write_reduced_all_csv(all_df, all_outpath):
    all_df = combined_resultFiles(abricate_inpath)

    to_drop = ["COVERAGE","COVERAGE_MAP","GAPS", "#FILE"]
    reduced_all_df = all_df.drop(columns=to_drop)
    reduced_all_df.rename(columns={"SEQUENCE": "CONTIG"}, inplace=True)
    reduced_all_df.to_csv(all_outpath)

path = "results/analysis/abricate/"
strain = "13"
abricate_inpath = "{}{}/*.tab".format(path, strain)
all_outpath = "{0}{1}/{1}_all.csv".format(path, strain)


#all_df = combined_resultFiles(abricate_inpath)
#write_reduced_all_csv(all_df, all_outpath)

all_df = pd.read_csv(all_outpath, index_col=0, header=0)

#sort by contig
per_contig_path = "{0}{1}/summarised/".format(path, strain) #need to create folder
contigs = all_df["CONTIG"].unique()
for contig in contigs:
    contig_df = all_df.loc[all_df["CONTIG"]==contig]
    sorted_contig_df = contig_df.sort_values(by=["START"])
    sorted_contig_df.to_csv("{0}/contig_{1}_sorted.csv".format(per_contig_path, contig))

#print(all_df.head)