#just once
import os
import glob

input_folder = "/homes/josefa/21-10-26_illumina/211026_VH00410_1_AAAGTLLHV/Analysis/3/Data/fastq/"
out_folder = "data/fastqs/"

i_files = glob.glob(f"{input_folder}JW*")

for i_file in i_files:
    os.system(f"gunzip -k {i_file}")
    unzip_i_file = i_file[:i_file.rfind(".gz")]
    print(unzip_i_file)
    os.system(f"mv {unzip_i_file} {out_folder}")

