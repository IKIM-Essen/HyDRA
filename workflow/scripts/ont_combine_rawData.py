## run once to have all ONT files correctly named in one folder
import glob
import os

infile = "/projects/seqlab/GridIONs/20211115_SM0006_JW0002_CefRes_strains_20211115/SM0006_JW0002_CefRes_strains_20211115/20211115_1648_X1_FAP91415_d0c3f758/fastq_pass/barcode10/FAP91415_pass_barcode10_f2d1671c_0.fastq.gz"
#infile = str(snakemake.input)
outfile = "data/ont/115.fastq.gz"
#outfile = str(snakemake.output)

search_in = infile[:infile.rfind("_")] + "*"
fastqs = glob.glob(search_in)

cmd = "cat"
for fastq in fastqs:
    cmd = "{0} {1}".format(cmd, fastq)
cmd = "{0} > {1}".format(cmd, outfile)

os.system(cmd)