import os
import glob
"""
strains = ["13", "62", "66", "102", "109", "115", "120", "129", "139", "156", "167", "169", "188", "265", "281"]

for strain in strains:
    nanoplot_path = f"results/reports/before_trim/nanoplot/{strain}/NanoPlot-report.html"
    new_np_path = f"results/reports/before_trim/nanoplot/{strain}_before_trim_NanoPlot-report.html"
    os.system(f"cp {new_np_path} {nanoplot_path}")

    nanoqc_path = f"results/reports/before_trim/nanoqc/{strain}/nanoQC.html"
    new_nq_path = f"results/reports/before_trim/nanoqc/{strain}_before_trim_nanoQC.html"
    os.system(f"cp {new_nq_path} {nanoqc_path}")
"""
#folders = glob.glob("old_filter_results/unicycler/*")
folders = glob.glob("results/unicycler/*")
download_folder = "to_download/"
for folder in folders:
    if folder.find("_") != -1:
        continue
    files = glob.glob(folder + "/*")
    strain = folder[folder.rfind("/") +1 :]
    for file in files:
        if file.find("assembly.gfa") != -1 and file.find("003") == -1:
            new_file = download_folder + strain + "_" + file[file.rfind("/") +1 :]
            os.system(f"cp {file} {new_file}")
        """
        if file.find("assembly.fasta") != -1:
            new_file = download_folder + strain + "_" + file[file.rfind("/") +1 :]
            os.system(f"cp {file} {new_file}")
        elif file.find("001_") != -1:
            new_file = download_folder + strain + "_" + file[file.rfind("/") +1 :]
            os.system(f"cp {file} {new_file}")
        """