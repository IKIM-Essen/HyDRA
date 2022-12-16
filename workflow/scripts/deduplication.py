"""
fastq_r1 = "dedup_test/13_R1.fastq"
fastq_r2 = "dedup_test/13_R2.fastq"

#fastq_r1 = "dedup_test/13_L001_R1_trimmed.fastq"
#fastq_r2 = "dedup_test/13_L001_R2_trimmed.fastq"

duplicates_path = fastq_r1[:fastq_r1.rfind(".")] + "_duplicates.txt"
dedupfastq1_path = fastq_r1[:fastq_r1.rfind(".")] + "_dedup.fastq"
dedupfastq2_path = fastq_r2[:fastq_r2.rfind(".")] + "_dedup.fastq"

"""
import os

def write_4_lines(file, header, seq, plus, quality):
    file.write(header)
    file.write(seq)
    file.write(plus)
    file.write(quality)

def deduplicate_r1(fastq_r1, duplicates_path, dedupfastq1_path):
    print("starting deduplication of R1")
    ffile = open(fastq_r1, "r")
    duplicates_file = open(duplicates_path, "w")
    newfastq = open(dedupfastq1_path, "w")

    num_duplicates =0
    starting_seqs = {}

    for line in ffile:
        if line.startswith("@VH"):
            header = line
            count = 2
            continue
        if count == 2:
            seq = line
            if seq in starting_seqs:
                duplicate = True
                num_duplicates += 1
            else:
                starting_seqs[seq] = header
                duplicate = False
            count += 1
            continue
        if duplicate:
            if count == 3:
                duplicates_file.write(header)
                count += 1
                continue
            else:
                continue
        else:
            if count == 3:
                plusline = line
                count += 1
                continue
            if count == 4:
                quality = line
                write_4_lines(newfastq, header, seq, plusline, quality)
    
    print(f"#duplicates: {num_duplicates}")
    print(f"#remaining: " + str(len(starting_seqs)))
  
    ffile.close()
    duplicates_file.close()
    newfastq.close()  

def create_duplicates_set(duplicates_path):
    duplicates_file = open(duplicates_path, "r")
    duplicates_set = set(())
    for line in duplicates_file:
        id = line.split()[0]
        duplicates_set.add(id)
    duplicates_file.close()
    return(duplicates_set)

def deduplicate_r2(dedupfastq2_path, fastq_r2, duplicates_path):
    print("starting deduplication of R2")
    newfastq_r2 = open(dedupfastq2_path, "w")
    f2file = open(fastq_r2, "r")
    duplicates_set = create_duplicates_set(duplicates_path)

    num_duplicates =0

    for line in f2file:
        if line.startswith("@VH"):
            id = line.split()[0]
            if id in duplicates_set:
                is_duplicate = True
                num_duplicates += 1
            else:
                is_duplicate = False
                newfastq_r2.write(line)
            continue
        if is_duplicate:
            continue
        else:
            newfastq_r2.write(line)
    
    print(f"#duplicates: {num_duplicates}")

    newfastq_r2.close()
    f2file.close()


fastq_r1 = str(snakemake.input.fastq1)
fastq_r2 = str(snakemake.input.fastq2)
dedupfastq1_path = str(snakemake.output.fastq1_dedup)
dedupfastq2_path = str(snakemake.output.fastq2_dedup)
duplicates_path = fastq_r1[:fastq_r1.rfind(".")] + "_duplicates.txt"

deduplicate_r1(fastq_r1, duplicates_path, dedupfastq1_path)
deduplicate_r2(dedupfastq2_path, fastq_r2, duplicates_path)
os.system(f"rm {duplicates_path}")