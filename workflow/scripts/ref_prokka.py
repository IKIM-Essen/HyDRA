import glob
import os
"""
refs = glob.glob("ressources/ref_genomes/GCF_ab*")
outdir_com = "results/analysis/ref_prokka/"
for ref in refs:
    prefix = ref[ref.rfind("/") +1 :ref.find("_0")]
    outdir = outdir_com + prefix
    os.system(f"prokka --outdir {outdir}/ --prefix {prefix} {ref}")
"""
out = "results/analysis/prokka/sum.txt"
sumfile = open(out, "w")
prokkas = glob.glob("results/analysis/prokka/*")

for prokka in prokkas:
    if prokka.find("mauve") >=0:
        continue
    elif prokka.find("sum") >=0:
        continue

    strain = prokka[prokka.rfind("/"):]
    sumline = ""
    txt_path = prokka + strain + ".txt"
    prokka_txt = open(txt_path, "r")
    for line in prokka_txt:
        #print(line)
        value = line[line.find(":") +2:]
        sumline = sumline + value.replace("\n", "\t")
    sumfile.write(prokka + "\n")
    sumfile.write(sumline + "\n")

    prokka_txt.close()
sumfile.close()
