import os
import glob

def minirmd(in_fastq1, in_fastq2, general_outfile, d_param, log):
    os.system(f"minirmd -i {in_fastq1} -f {in_fastq2} -o {general_outfile} -d {d_param} > {log} 2>&1")

def renaming(general_outfile, out_fastq1, out_fastq2):
    rmd_outfiles = glob.glob(general_outfile + "*")
    for rmd_outfile in rmd_outfiles:
        suffix = rmd_outfile[rmd_outfile.rfind("_"):]
        if suffix == "_1":
            os.system(f"mv {rmd_outfile} {out_fastq1}")
        elif suffix == "_2":
            os.system(f"mv {rmd_outfile} {out_fastq2}")

in_fastq1 = snakemake.input[0]
in_fastq2 = snakemake.input[1]
d_param = snakemake.params
log = snakemake.log

out_fastq1 = snakemake.output[0]
out_fastq2 = snakemake.output[1]

general_outfile = out_fastq1[:out_fastq1.rfind("_R")]

minirmd(in_fastq1, in_fastq2, general_outfile, d_param, log)
renaming(general_outfile, out_fastq1, out_fastq2)
