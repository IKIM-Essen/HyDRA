import argparse
from os import listdir
import pandas as pd

ASSEMBLY_TYPES = ["sample_name", "long", "shortR1" ,"shortR2" ,"assembly"]

df_samples = pd.DataFrame(columns=ASSEMBLY_TYPES)

def generate_sample_file(assembly_type, path, output):
    index = 0
    for file_name in listdir(path):
        df_samples.loc[index, "sample_name"] = file_name.replace(".fna.gz", "")
        if(assembly_type == "none"):
            df_samples.loc[index, "assembly"] = path + "/" + file_name
        elif(assembly_type == "short" or assembly_type == "hybrid"):
            file_name_r1 = file_name.replace(".fna.gz", "_R1.fna.gz")
            df_samples.loc[index, "shortR1"] = path + "/" + file_name_r1
            file_name_r2 = file_name.replace(".fna.gz", "_R2.fna.gz")
            df_samples.loc[index, "shortR2"] = path + "/" + file_name_r2
        if(assembly_type == "long" or assembly_type == "hybrid"):
            df_samples.loc[index, "long"] = path + "/" + file_name


        index = index + 1
        
    df_samples.to_csv(output, sep=',', index=False)

parser = argparse.ArgumentParser(description="Generate samples.csv for all files in folder")
parser.add_argument("folder_input", help="Path to the input folder")
parser.add_argument("assembly_type_input", help="Choose one of these assembly types: " + "short | long | hybrid | none")
parser.add_argument("file_output", help="Output file")
args = parser.parse_args()
generate_sample_file(args.assembly_type_input, args.folder_input, args.file_output)
print(args.file_output + " created!")