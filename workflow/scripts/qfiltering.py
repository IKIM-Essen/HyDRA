"""
quality = 26 #snakemake.params.q
percentage = 80 #snakemake.params.p
inputfile1 = "data/illumina/13_L002_R1.fastq" #snakemake.input.fastq1
inputfile2 = "data/illumina/13_L002_R2.fastq" #snakemake.input.fastq2
filtered1 = "results/preprocess_ill/13/13_L002_R1_filtered.fastq" #snakemake.output.fastq1_filt
filtered2 = "results/preprocess_ill/13/13_L002_R2_filtered.fastq" #snakemake.output.fastq2_filt
unpairedfile = "results/preprocess_ill/13/13_L002_unpaired.fastq" #snakemake.output.unpaired
"""
import os

def qfiltering(quality, percentage, inputfile, outputfile, logfile):
    os.system(f"fastq_quality_filter -q {quality} -p {percentage} -i {inputfile} -o {outputfile} > {logfile} 2>&1")

def make_header_set(inputfile):
    ffile = open(inputfile, "r")
    header_set = set(())
    for line in ffile:
        if line.startswith("@VH"):
            header_set.add(line)
        else:
            continue
    ffile.close()
    return(header_set)

def compare_sets(header_unfiltered, header_filtered):
    outfiltered_set = set(())
    for header in header_unfiltered:
        if header not in header_filtered:
            id = header.split()[0]
            outfiltered_set.add(id)
    return(outfiltered_set)

def missing_in_one(outfiltered1, outfiltered2):
    missing_dict = {}
    for header_id in outfiltered1:
        if header_id not in outfiltered2:
            missing_dict[header_id] = 1
    for header_id in outfiltered2:
        if header_id not in outfiltered1:
            missing_dict[header_id] = 2
    return(missing_dict)

def write_outfiles(prefiltered1, prefiltered2, filtered1, filtered2, unpairedfile, missing_in_one_dict):
    out1 = open(filtered1, "w")
    out2 = open(filtered2, "w")
    unpaired = open(unpairedfile, "w")
    pre1 = open(prefiltered1, "r")
    pre2 = open(prefiltered2, "r")

    for line in pre1:
        if line.startswith("@VH"):
            id = line.split()[0]
            if id in missing_in_one_dict:
                if missing_in_one_dict[id] == 2:
                    is_unpaired = True
                    unpaired.write(line)
            else:
                is_unpaired = False
                out1.write(line)
        else:
            if is_unpaired:
                unpaired.write(line)
            else:
                out1.write(line)
    
    for line in pre2:
        if line.startswith("@VH"):
            id = line.split()[0]
            if id in missing_in_one_dict:
                if missing_in_one_dict[id] == 1:
                    is_unpaired = True
                    unpaired.write(line)
            else:
                is_unpaired = False
                out2.write(line)
        else:
            if is_unpaired:
                unpaired.write(line)
            else:
                out2.write(line)            

    pre1.close()
    pre2.close()
    out1.close()
    out2.close()
    unpaired.close()


quality = str(snakemake.params.q)
percentage = str(snakemake.params.p)
inputfile1 = str(snakemake.input.fastq1)
inputfile2 = str(snakemake.input.fastq2)
filtered1 = str(snakemake.output.fastq1_filt)
filtered2 = str(snakemake.output.fastq2_filt)
unpairedfile = str(snakemake.output.unpaired)
logfile = str(snakemake.log)
prefiltered1 = filtered1.replace("filtered", "prefiltered")
prefiltered2 = filtered2.replace("filtered", "prefiltered")


header_unfiltered1 = make_header_set(inputfile1)
print("#unfiltered sequences R1: " + str(len(header_unfiltered1)))
header_unfiltered2 = make_header_set(inputfile2)
print("#unfiltered sequences R2: " + str(len(header_unfiltered2)))

qfiltering(quality, percentage, inputfile1, prefiltered1, logfile)
qfiltering(quality, percentage, inputfile2, prefiltered2, logfile)

header_filtered1 = make_header_set(prefiltered1)
print("#filtered sequences R1: " + str(len(header_filtered1)))
header_filtered2 = make_header_set(prefiltered2)
print("#filtered sequences R2: " + str(len(header_filtered2)))

outfiltered1 = compare_sets(header_unfiltered1, header_filtered1)
print("#sequences R1 filtered out: " + str(len(outfiltered1)))
outfiltered2 = compare_sets(header_unfiltered2, header_filtered2)
print("#sequences R2 filtered out: " + str(len(outfiltered2)))

missing_in_one_dict = missing_in_one(outfiltered1, outfiltered2)
print("#unpaired sequences: " + str(len(missing_in_one_dict)))

write_outfiles(prefiltered1, prefiltered2, filtered1, filtered2, unpairedfile, missing_in_one_dict)
os.system(f"rm {prefiltered1} {prefiltered2}")