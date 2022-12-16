import pandas as pd
import numpy as np

def calculate_gc_content(seq):
    g_and_c = seq.count("G") + seq.count("C")
    content = g_and_c / len(seq) * 100
    return(round(content, 2))

def calculate_gene_gc_skew(seq):
    g_count = seq.count("G")
    c_count = seq.count("C")
    skew = (g_count - c_count)/(g_count + c_count)
    return(round(skew, 2))

def create_allGenes_file(gff_path):
    gfile = open(gff_path + ".gff", "r")
    all_genes_file = open(gff_path + "_allGenes.gff", "w")

    for line in gfile:
        if line.startswith("##FASTA"):
            break
        if line.startswith("##"):
            continue
        
        all_genes_file.write(line)

    gfile.close()
    all_genes_file.close()

def number_identifier_properties(ident_list):
    properties_list = []
    for list in ident_list:
        for entry in list:
            property = entry.split("=")[0]
            if property not in properties_list:
                properties_list.append(property)

    print(properties_list)
    print("number of different values: " + str(len(properties_list)))

def check_property_equality(ident_list):
    for list in ident_list:
        for entry in list:    
            property = entry.split("=")
            if property[0] == "ID":
                id = property[-1]
                continue
            if property[0] == "locus_tag":
                if id != property[-1]:
                    print("ID and locus_tag not the same!")
                continue
            if property[0] == "Name":
                genename = property[-1]
                continue
            if property[0] == "gene":
                if genename != property[-1]:
                    print("Name and gene not the same!")
                continue

def reduce_inference(inference):
    comma_pos = inference.find(",")
    if comma_pos >= 0:
        newinf = inference[comma_pos +1 :]
        return(newinf[newinf.find(":") +1 :])
    else:
        return("")

def seperate_info_string(gff_path):
    all_genes_file = open(gff_path + "_allGenes.gff", "r")
    all_genes_tab = open(gff_path + "_allGenes.tab", "w")
    ident_list = []
    properties_list = ["ID", "gene", "product", "db_xref", "eC_number", "inference", "note"]

    header = "contig\ttool\tftype\tstart\tend\tp-value\tstrand\tframe\t" + "\t".join(properties_list)
    all_genes_tab.write(header)
    for line in all_genes_file:
        split_line = line.split("\t")
        identifier = split_line[-1]

        split_ident = identifier.split(";")
        ident_list.append(split_ident)

        properties_dict = {}
        value_list = []

        for property in split_ident:
            properties_dict[property.split("=")[0]] = property.split("=")[1].strip()
        for property in properties_list:
            if property in properties_dict:
                if property == "inference":
                    value_list.append(reduce_inference(properties_dict[property]))
                elif property == "db_xref":
                    colon_index = properties_dict[property].find(":")
                    value_list.append(properties_dict[property][colon_index +1 :])
                else:
                    value_list.append(properties_dict[property])
            else:
                value_list.append("")

        split_line.pop(-1)
        split_line.extend(value_list)
        new_line = "\t".join(split_line)
        all_genes_tab.write("\n" + new_line)
    all_genes_file.close()
    all_genes_tab.close()

def write_out(df, csv_path, tab_path):
    #writes csv
    df.to_csv(path_or_buf=csv_path, float_format="%.2f")
    #writes tab seperated file
    df.to_csv(path_or_buf=tab_path, sep='\t', float_format="%.2f")

def create_gc_df(ffn_path):
    ffn_file = open(gff_path + ".ffn", "r")
    sequence_dict = {}
    sequence = ""
    for line in ffn_file:
        if line.startswith(">"):
            if sequence != "":
                sequence_dict[id] = sequence
            id = line.split()[0][1:]
            sequence = ""
            continue
        else:
            sequence += line.strip()
    sequence_dict[id] = sequence
    ffn_file.close()

    gc_dict = {}
    #calculation of gc content & skew
    for id in sequence_dict:
        gc_content = calculate_gc_content(sequence_dict[id])
        gc_skew = calculate_gene_gc_skew(sequence_dict[id])
        gc_dict[id] = [gc_content, gc_skew]

    gc_df = pd.DataFrame.from_dict(gc_dict, orient='index', columns=["gc_content", "gc_skew"])
    return(gc_df)

def add_gc_values(ffn_path, test_df):
    #add gc content & skew to dataframe
    gc_df = create_gc_df(ffn_path)

    final_df = test_df.join(gc_df)
    final_df.insert(8,"gc_content", final_df.pop("gc_content"))
    final_df.insert(9,"gc_skew", final_df.pop("gc_skew"))
    return(final_df)

def get_value_overview(column_list):
    for column in column_inspection:
        dif_values = test_df[column].unique()
        print("There are {0} different {1} values.\t{2}".format(str(len(dif_values)), column, str(dif_values)))


gff_path = "ml_ABR/115"
#create_allGenes_file(gff_path)
#seperate_info_string(gff_path)

test_df = pd.read_table(gff_path + "_allGenes.tab", header=0, index_col="ID")

column_inspection = ["contig", "ftype", "strand", "frame"]
get_value_overview(column_inspection)

ffn_path = gff_path + ".ffn"

gc_df = add_gc_values(ffn_path, test_df)

#length calculation
gc_df.insert(3, "length", test_df["end"] - test_df["start"]+1)

#drop redundant & unimportant features
list_to_reduce = ["tool", "frame", "p-value", "start", "end", "note"]
reduced_df = gc_df.drop(columns=list_to_reduce)

#convert ftype & strand into numeric
ftype_dict = {'CDS':0, 'rRNA':1, 'tRNA':2, 'tmRNA':3}
strand_dict = {'+':0, '-':1}
reduced_df = reduced_df.replace({"ftype":ftype_dict}).replace({"strand":strand_dict})


output_path = gff_path[:gff_path.find("/")] + "/results" + gff_path[gff_path.find("/"):]
#writes all genes to files (tab delemited and csv)
csv_path = output_path + "_allcontigs.csv"
tab_path = output_path + "_allcontigs.tab"
write_out(reduced_df, csv_path, tab_path)

#writes all genes for one contig to file (tab delemited and csv)
contigs = reduced_df["contig"].unique()
for contig in contigs:
    contig_df = reduced_df.loc[reduced_df["contig"] == contig]
    contig_df = contig_df.drop(columns="contig")
    csv_path = output_path + "_contig{0}.csv".format(contig)
    tab_path = output_path + "_contig{0}.tab".format(contig)
    write_out(contig_df, csv_path, tab_path)


"""
cogs = reduced_df["inference"].unique()
for cog in cogs:
    if type(cog)!= str:
        if np.isnan(cog):
            continue
    if cog.count(":")>=2:
        continue
    else:
        print(cog)
        break
print("done")


all_genes_tab = open(gff_path + "_allGenes.tab", "r")

ident_list = []
inf_list = []
for line in all_genes_tab:
    split_line = line.split("\t")
    
    #remove unnecessary information
    inference = split_line[-2]
    comma_pos = inference.find(",")
    if comma_pos >= 0:
        split_line[-2] = inference[comma_pos +1 :]
    else:
        split_line[-2] = ""
    
    end = 
    length = end - start +1
    
    split_ident = identifier.split(";")
    ident_list.append(split_ident)

    properties_dict = {}
    value_list = []

    for property in split_ident:
        properties_dict[property.split("=")[0]] = property.split("=")[1].strip()
    for property in properties_list:
        if property in properties_dict:
            value_list.append(properties_dict[property])
        else:
            value_list.append("")

    split_line.pop(-1)
    split_line.extend(value_list)
    new_line = "\t".join(split_line)
    all_genes_tab.write(new_line + "\n")


    #tool = split_line[1]
 
print(inf_list)
all_genes_tab.close()
#number_identifier_properties(ident_list)
#check_property_equality(ident_list)
"""