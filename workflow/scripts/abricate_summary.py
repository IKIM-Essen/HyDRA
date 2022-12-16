import glob
import os

def abricate(db, assembly, strain, outdir, logfile):
    os.system(f"abricate --db {db} {assembly} > {outdir}/{strain}_{db}.tab 2>> {logfile}") #funktionierts mit 2 >?

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

dbs = ["ncbi", "plasmidfinder", "argannot", "megares", "ecoli_vf", "resfinder", "card", "ecoh", "vfdb"]

assembly = str(snakemake.input) #"results/final_assemblies/115_copy/assembly.fasta" #str(snakemake.input)
outdir = str(snakemake.params.outdir) #"results/analysis/abricate/115/" #str(snakemake.params.outdir)
strain = str(snakemake.wildcards.strain)
logfile = str(snakemake.log)

for db in dbs:
    abricate(db, assembly, strain, outdir, logfile)
write_summary(outdir, strain)
