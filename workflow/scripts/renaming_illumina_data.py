import glob
import os

to_rename_dict = {"1":"91","2":"102","3":"109","4":"129","5":"139",
    "6":"156","7":"167","8":"169","9":"188","vs17687":"281"}
to_rename = to_rename_dict.keys()

tfile = open("data/renaming_illumina.txt","a")
search_path = "data/fastqs/"

fastqs = glob.glob(f"{search_path}*")
for fastq in fastqs:

    fname = fastq.split("/")[-1]

    details = fname.split("_")
    strain = details[1]
    lane = details[3]
    read = details[4]

    if strain in to_rename:
        strain = to_rename_dict[strain]

    new_name = f"{strain}_{lane}_{read}.fastq"
    new_path = search_path + new_name

    #write renaming in file
    tfile.write(fname + "\n")
    tfile.write(new_name + "\n\n")

    os.system(f"mv {fastq} {new_path}")

tfile.close()