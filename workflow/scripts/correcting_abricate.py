import glob
import pandas as pd

#combines the results of the files (eg 115_ncbi.tab) in path folder
#drops file & coverage map columns + renames column sequence in contig
#saves all results in one tab file at all_path
def combine_results(path, all_path):

    tab_paths = glob.glob(path + "115*.tab")
    all_file = open(all_path, "w")

    first_tabfile = True

    for tab_path in tab_paths:
        tab_file = open(tab_path, "r")
        for line in tab_file:
            sp_line = line.split("\t")
            if first_tabfile:
                index_covmap = sp_line.index("COVERAGE_MAP")
                sp_line[sp_line.index("SEQUENCE")] = "CONTIG"
                first_tabfile = False
            elif line.startswith("#"):
                continue
            
            sp_line.pop(index_covmap)
            sp_line.pop(0)

            all_file.write("\t".join(sp_line))
        
        tab_file.close()
    all_file.close()

path = "results/analysis/abricate/115/"
all_path = path + "all_115.tab"
#combine_results(path, all_path)

test_df = pd.read_table(all_path, header=0)
sorted = test_df.sort_values(by="CONTIG")

contigs_list = []
for entry in sorted['CONTIG']:
    if entry not in contigs_list:
        contigs_list.append(entry)

sorted_csv_path = path + "115_contig__sorted.csv"
insertion_index = sorted_csv_path.find("_sort")
for contig in contigs_list:
    filtered_for_contig = sorted.loc[lambda df: df['CONTIG'] == contig, :]
    sorted_contig = filtered_for_contig.sort_values(by="START")
    sorted_contig_csv_path = sorted_csv_path[: insertion_index] + str(contig) + sorted_csv_path[insertion_index :]
    sorted_contig.to_csv(path_or_buf=sorted_contig_csv_path, index=False)


"""

#sorted_contigs = sorted['CONTIG']

sorted_contig2 = contig2.sort_values(by="START")
#sorted_contig2.to_csv(path_or_buf=sorted_csv_path)
#print(sorted_contig2)

#print(sorted_contigs.loc[14:24])
#print(sorted_contigs)

#print(sorted.loc['1','CONTIG'])
#sorted.to_csv(path_or_buf=sorted_csv_path)


#sort by contig
all_file = open(all_path, "r")
first_line = True
lines_list = []
for line in all_file:
    if first_line:
        first_line = False
        continue
    sp_line = line.split("\t")
    lines_list.append(sp_line)
all_file.close()

for entry in lines_list:
    
  """  

"""
sumfile_path = path + "115_summary.txt"
sumfile = open(sumfile_path, "r")

for line in sumfile:
    if line.startswith("genes"):
    
        print(line.split()[-1])

"""