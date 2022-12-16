import os
import glob

input_folder = "results/flye/{id}_{stage}/" #from snakemake flye output
output_folder = "results/checkm/{id}_{stage}/" #from snakemake checkm output
marker_file = "results/checkm/marker_file"

checkm_path = "results/reports/assembly/checkm/"
checkm_folders = glob.glob(checkm_path + "*")
for folder in checkm_folders:
    strain = folder[folder.rfind("/") +1:]
    log_file = f"logs/checkm/final_{strain}.log"
    os.system(f"cp {log_file} {folder}")