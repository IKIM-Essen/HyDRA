import os

#filters all genes with hypthetical protein as product/description out
def filter_hypothetical_proteins(input_tsv):
    number_hypo = 0
    wo_hyprot_path = input_tsv.replace(".tsv", "_wo_hyprot.tsv")
    if wo_hyprot_path.find("ref_prokka") >= 0:
        wo_hyprot_path = wo_hyprot_path.replace("ref_prokka", "ref_vs_prokka")
    else:
        wo_hyprot_path = wo_hyprot_path.replace("prokka", "ref_vs_prokka")
    tsv = open(input_tsv, "r")
    wo_hyprot = open(wo_hyprot_path, "w")
    for gene in tsv:
        description = gene.split("\t")[-1].strip()
        if description != "hypothetical protein":
            wo_hyprot.write(gene)
        else:
            number_hypo += 1

    tsv.close()
    wo_hyprot.close()
    print("#hypothetical proteins: " + str(number_hypo))

#seperates the .tsv file by ftype into cds, cds without genename, rna
def seperate_ftypes(outpath, ref_tsv_path, ref_name):
    number_cds = number_tmrna = number_trna = number_rrna = number_cds_wo_gname = 0
    tsv = open(ref_tsv_path, "r")

    cds_path = f"{outpath}{ref_name}/{ref_name}_cds.tsv"
    cds_file = open(cds_path, "w")
    cds_wo_gname_path = f"{outpath}{ref_name}/{ref_name}_cds_wog.tsv"
    cds_wo_gname_file = open(cds_wo_gname_path, "w")
    rna_path = f"{outpath}{ref_name}/{ref_name}_rna.tsv"
    rna_file = open(rna_path, "w")

    for gene in tsv:
        splitted = gene.split("\t")
        ftype = splitted[1]
        is_rna = False
        if ftype == "CDS":
            gene_name = splitted[3]
            if gene_name == "":
                number_cds_wo_gname += 1
                cds_wo_gname_file.write(gene)
            else:
                number_cds += 1
                cds_file.write(gene)
        elif ftype == "tmRNA":
            number_tmrna += 1
            is_rna = True
        elif ftype == "rRNA":
            number_rrna += 1
            is_rna = True
        elif ftype == "tRNA":
            number_trna += 1
            is_rna = True
        if is_rna:
            rna_file.write(gene)

    cds_file.close()
    rna_file.close()
    cds_wo_gname_file.close()

    print(f"{ref_name}:")
    print(f"#CDS without genename: {number_cds_wo_gname}")
    print(f"#CDS: {number_cds}")
    print(f"#tmRNA: {number_tmrna}")
    print(f"#tRNA: {number_trna}")
    print(f"#rRNA: {number_rrna}")

#returns a list of dictionaries for the different types of rna (dict{description: gene_line from tsv file}) 
def make_rna_dicts(rna_tsv_path):
    trna_dict = {}
    rrna_dict = {}
    tmrna_dict = {}

    rna_tsv_file = open(rna_tsv_path, "r")
    for gene in rna_tsv_file:
        splitted = gene.split("\t")
        ftype = splitted[1]
        gene_name = splitted[3]
        description = splitted[-1].strip()
        if ftype == "tmRNA":
            tmrna_dict[gene_name] = gene
        elif ftype == "tRNA":
            if description in trna_dict:
                trna_dict[description].append(gene)
            else:
                trna_dict[description] = [gene]
        elif ftype == "rRNA":
            if description in rrna_dict:
                rrna_dict[description].append(gene)
            else:
                rrna_dict[description] = [gene]
    rna_tsv_file.close()

    print("#different tmRNA: " + str(len(tmrna_dict)))
    print("#different tRNA: " + str(len(trna_dict)))
    print("#different rRNA: " + str(len(rrna_dict)))

    return([tmrna_dict, trna_dict, rrna_dict])

#creates dictionaries for genes with genename as key
def make_cds_dict(cds_tsv_path):
    cds_dict = {}

    cds_tsv_file = open(cds_tsv_path, "r")
    for gene in cds_tsv_file:
        splitted = gene.split("\t")
        ftype = splitted[1]
        gene_name = splitted[3]
        description = splitted[-1].strip()
        if ftype == "CDS":
            if gene_name in cds_dict:
                cds_dict[gene_name].append(gene)
            else:
                cds_dict[gene_name] = [gene]
    cds_tsv_file.close()

    print("#different CDS in reference: " + str(len(cds_dict)))
    return(cds_dict)

#compares dictionaries of strain & reference of genes and writes output files
def cds_ref_vs_strain(cds_dict, strain_cds_dict, outpath_strain):
    cds_both_dict = {"different":[], "same":[]}
    cds_ref_lonly_list = []
    for ref_genename in cds_dict:
        if ref_genename in strain_cds_dict:
            refgene_list = cds_dict[ref_genename]
            straingene_list = strain_cds_dict[ref_genename]
            #easiest case
            if len(refgene_list) == 1 and len(straingene_list) == 1:
                refg_linesplit = refgene_list[0].split("\t")
                straing_linesplit = straingene_list[0].split("\t")
                refg_len = refg_linesplit[2]
                straing_len = straing_linesplit[2]
                if refg_len == straing_len:
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(refg_linesplit[3], refg_linesplit[-1].strip(), refg_len, refg_linesplit[0], straing_linesplit[0])
                    cds_both_dict["same"].append(outstr)
                else:
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(refg_linesplit[3], refg_linesplit[-1].strip(), refg_len, straing_len, refg_linesplit[0], straing_linesplit[0])
                    cds_both_dict["different"].append(outstr)
            else:
                print("take closer look! gene lists longer than 1")
        else:
            if len(cds_dict[ref_genename]) > 1:
                print("take a look! only in ref list longer 1")
            cds_ref_lonly_list.append(cds_dict[ref_genename][0])

    cds_strain_lonly_list = []
    for ref_genename in strain_cds_dict:
        if ref_genename not in cds_dict:
            if len(strain_cds_dict[ref_genename]) > 1:
                print("take a look! only in strain list longer 1")
            cds_strain_lonly_list.append(strain_cds_dict[ref_genename][0])

    print("#same genes: " + str(len(cds_both_dict["same"])))
    cds_both_same_path = f"{outpath_strain}_cds_both_same.tsv"
    header = "gene\tproduct\tref_length_bp\tref_locus_tag\tstrain_locus_tag\n"
    write_out_comparison(cds_both_same_path, cds_both_dict["same"], header)

    print("#different genes: " + str(len(cds_both_dict["different"])))
    cds_both_diff_path = f"{outpath_strain}_cds_both_diff.tsv"
    header = "gene\tproduct\tref_length_bp\tstrain_length_bp\tref_locus_tag\tstrain_locus_tag\n"
    write_out_comparison(cds_both_diff_path, cds_both_dict["different"], header)

    print("#genes only in strain: " + str(len(cds_strain_lonly_list)))
    cds_only_strain_path = f"{outpath_strain}_cds_only_strain.tsv"
    header = "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG	product\n"
    write_out_comparison(cds_only_strain_path, cds_strain_lonly_list, header)

    print("#genes only in reference: " + str(len(cds_ref_lonly_list)))
    cds_only_ref_path = f"{outpath_strain}_cds_only_ref.tsv"
    header = "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG	product\n"
    write_out_comparison(cds_only_ref_path, cds_ref_lonly_list, header)

def write_out_comparison(file_path, outstr_list, header):
    ffile = open(file_path, "w")
    ffile.write(header)
    for entry in outstr_list:
        ffile.write(entry)
    ffile.close()

#creates dictionaries for genes without genenames, with COG, EC number or description as key
def make_cds_wog_dicts(cds_wog_tsv_path):
    cds_wog_cogs_dict = {}
    cds_wog_ec_dict = {}
    cds_wog_only_des_dict = {}

    cds_wog_tsv_file = open(cds_wog_tsv_path, "r")
    number_cogs = number_ec = number_des = 0
    for gene in cds_wog_tsv_file:
        splitted = gene.split("\t")
        if splitted[-2] != "":
            cog = splitted[-2]
            number_cogs += 1
            if cog in cds_wog_cogs_dict:
                cds_wog_cogs_dict[cog].append(gene)
            else:
                cds_wog_cogs_dict[cog] = [gene]
        elif splitted[-3] != "":
            ec = splitted[-3]
            number_ec += 1
            if ec in cds_wog_ec_dict:
                cds_wog_ec_dict[ec].append(gene)
            else:
                cds_wog_ec_dict[ec] = [gene]
        else:
            des = splitted[-1]
            number_des += 1
            if des in cds_wog_ec_dict:
                cds_wog_only_des_dict[des].append(gene)
            else:
                cds_wog_only_des_dict[des] = [gene]

    print(f"#COG genes {number_cogs}")
    print(f"#EC genes {number_ec}")
    print(f"#genes without genename, COG & EC: {number_des}")
    cds_wog_tsv_file.close()
    return([cds_wog_cogs_dict, cds_wog_ec_dict, cds_wog_only_des_dict])

#compares dictionaries of strain & reference of genes without genename and writes output files
def cds_wog_ref_vs_strain(cds_wog_dict_list, cds_wog_strain_dict_list, outpath_strain):
    cds_wog_both_dict = {"different":[], "same":[]}
    cds_wog_ref_only_list = []
    cds_wog_strain_only_list = []
    uncommon_list = []
    for ref_cog in cds_wog_dict_list[0]:
        if ref_cog in cds_wog_strain_dict_list[0]:
            refgene_list = cds_wog_dict_list[0][ref_cog]
            straingene_list = cds_wog_strain_dict_list[0][ref_cog]
            #easiest case
            if len(refgene_list) == 1 and len(straingene_list) == 1:
                refg_linesplit = refgene_list[0].split("\t")
                straing_linesplit = straingene_list[0].split("\t")
                refg_len = refg_linesplit[2]
                straing_len = straing_linesplit[2]
                if refg_len == straing_len:
                    #EC COG description length id_ref id_strain
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(refg_linesplit[-3], refg_linesplit[-2], refg_linesplit[-1].strip(), refg_len, refg_linesplit[0], straing_linesplit[0])
                    cds_wog_both_dict["same"].append(outstr)
                else:
                    #EC COG description length_ref length_strain id_ref id_strain
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(refg_linesplit[-3], refg_linesplit[-2], refg_linesplit[-1].strip(), refg_len, straing_len, refg_linesplit[0], straing_linesplit[0])
                    cds_wog_both_dict["different"].append(outstr)
            else:
                for refgene in refgene_list:
                    uncommon_list.append(refgene)
                for straingene in straingene_list:
                    uncommon_list.append(straingene)
        else:
            if len(cds_wog_dict_list[0][ref_cog]) > 1:
                for element in cds_wog_dict_list[0][ref_cog]:
                    cds_wog_ref_only_list.append(element)
                #print("take a look! only in ref wog_cog_list longer 1")
            cds_wog_ref_only_list.append(cds_wog_dict_list[0][ref_cog][0])

    for strain_cog in cds_wog_strain_dict_list[0]:
        if strain_cog not in cds_wog_dict_list[0]:
            if len(cds_wog_strain_dict_list[0][strain_cog]) > 1:
                for element in cds_wog_strain_dict_list[0][strain_cog]:
                    cds_wog_strain_only_list.append(element)
                    #print("take a look! only in strain wog_cog_list longer 1")
            cds_wog_strain_only_list.append(cds_wog_strain_dict_list[0][strain_cog][0])

    #EC
    for ref_ec in cds_wog_dict_list[1]:
        if ref_ec in cds_wog_strain_dict_list[1]:
            refgene_list = cds_wog_dict_list[1][ref_ec]
            straingene_list = cds_wog_strain_dict_list[1][ref_ec]
            #easiest case
            if len(refgene_list) == 1 and len(straingene_list) == 1:
                refg_linesplit = refgene_list[0].split("\t")
                straing_linesplit = straingene_list[0].split("\t")
                refg_len = refg_linesplit[2]
                straing_len = straing_linesplit[2]
                if refg_len == straing_len:
                    #EC COG description length id_ref id_strain
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(refg_linesplit[-3], refg_linesplit[-1].strip(), refg_len, refg_linesplit[0], straing_linesplit[0])
                    cds_wog_both_dict["same"].append(outstr)
                else:
                    #EC COG description length_ref length_strain id_ref id_strain
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(refg_linesplit[-3], refg_linesplit[-1].strip(), refg_len, straing_len, refg_linesplit[0], straing_linesplit[0])
                    cds_wog_both_dict["different"].append(outstr)
            else:
                for refgene in refgene_list:
                    uncommon_list.append(refgene)
                for straingene in straingene_list:
                    uncommon_list.append(straingene)
        else:
            if len(cds_wog_dict_list[1][ref_ec]) > 1:
                print("take a look! only in ref wog_cog_list longer 1")
            cds_wog_ref_only_list.append(cds_wog_dict_list[1][ref_ec][0])

    for strain_ec in cds_wog_strain_dict_list[1]:
        if strain_ec not in cds_wog_dict_list[1]:
            if len(cds_wog_strain_dict_list[1][strain_ec]) > 1:
                    print("take a look! only in strain wog_ec_list longer 1")
            cds_wog_strain_only_list.append(cds_wog_strain_dict_list[1][strain_ec][0])

    #description
    for ref_des in cds_wog_dict_list[2]:
        if ref_des in cds_wog_strain_dict_list[2]:
            refgene_list = cds_wog_dict_list[2][ref_des]
            straingene_list = cds_wog_strain_dict_list[2][ref_des]
            #easiest case
            if len(refgene_list) == 1 and len(straingene_list) == 1:
                refg_linesplit = refgene_list[0].split("\t")
                straing_linesplit = straingene_list[0].split("\t")
                refg_len = refg_linesplit[2]
                straing_len = straing_linesplit[2]
                if refg_len == straing_len:
                    #EC COG description length id_ref id_strain
                    outstr = "{0}\t{1}\t{2}\t{3}\n".format(refg_linesplit[-1].strip(), refg_len, refg_linesplit[0], straing_linesplit[0])
                    cds_wog_both_dict["same"].append(outstr)
                else:
                    #EC COG description length_ref length_strain id_ref id_strain
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(refg_linesplit[-1].strip(), refg_len, straing_len, refg_linesplit[0], straing_linesplit[0])
                    cds_wog_both_dict["different"].append(outstr)
            else:
                for refgene in refgene_list:
                    uncommon_list.append(refgene)
                for straingene in straingene_list:
                    uncommon_list.append(straingene)
        else:
            if len(cds_wog_dict_list[2][ref_des]) > 1:
                print("take a look! only in ref wog_cog_list longer 1")
            cds_wog_ref_only_list.append(cds_wog_dict_list[2][ref_des][0])

    for strain_des in cds_wog_strain_dict_list[2]:
        if strain_des not in cds_wog_dict_list[2]:
            if len(cds_wog_strain_dict_list[2][strain_des]) > 1:
                    print("take a look! only in strain wog_des_list longer 1")
            cds_wog_strain_only_list.append(cds_wog_strain_dict_list[2][strain_des][0])

    print("#same wog_genes: " + str(len(cds_wog_both_dict["same"])))
    cds_wog_both_same_path = f"{outpath_strain}_cds_wog_both_same.tsv"
    header = "EC\tCOG\tproduct\tref_length_bp\tref_locus_tag\tstrain_locus_tag\n"
    write_out_comparison(cds_wog_both_same_path, cds_wog_both_dict["same"], header)

    print("#genes only in strain: " + str(len(cds_wog_strain_only_list)))
    cds_wog_only_strain_path = f"{outpath_strain}_cds_wog_only_strain.tsv"
    header = "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG	product\n"
    write_out_comparison(cds_wog_only_strain_path, cds_wog_strain_only_list, header)

    print("#genes only in reference: " + str(len(cds_wog_ref_only_list)))
    cds_wog_only_ref_path = f"{outpath_strain}_cds_wog_only_ref.tsv"
    header = "locus_tag\tftype\tlength_bp\tgene\tEC_number\tCOG	product\n"
    write_out_comparison(cds_wog_only_ref_path, cds_wog_ref_only_list, header)

    print("#different genes: " + str(len(cds_wog_both_dict["different"])))
    print("#uncommon genes: " + str(len(uncommon_list)))
    cds_wog_both_diff_path = f"{outpath_strain}_cds_wog_both_diff.tsv"
    cds_wog_both_diff_file = open(cds_wog_both_diff_path, "w")
    for entry in cds_wog_both_dict["different"]:
        cds_wog_both_diff_file.write(entry)
    cds_wog_both_diff_file.write("\n")
    for entry in uncommon_list:
        cds_wog_both_diff_file.write(entry)
    cds_wog_both_diff_file.close()

#compares RNA dictionaries of strain & reference and writes output files
def rna_ref_vs_strain(rna_dicts_list, rna_strain_dicts_list, outpath_strain):
    rna_both_dict = {"different":[], "same":[]}
    rna_ref_only_list = []
    rna_strain_only_list = []
    uncommon_list = []
    #tmRNA
    for ref_genename in rna_dicts_list[0]:
        if ref_genename in rna_strain_dicts_list[0]:
            refg_linesplit = rna_dicts_list[0][ref_genename].split("\t")
            straing_linesplit = rna_strain_dicts_list[0][ref_genename].split("\t")
            refg_len = refg_linesplit[2]
            straing_len = straing_linesplit[2]
            if refg_len == straing_len:
                #genename description length_ref id_ref id_strain
                outstr = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(refg_linesplit[3], refg_linesplit[-1].strip(), refg_len, refg_linesplit[0], straing_linesplit[0])
                rna_both_dict["same"].append(outstr)
            else:
                #genename description length_ref length_strain id_ref id_strain
                outstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(refg_linesplit[3], refg_linesplit[-1].strip(), refg_len, straing_len, refg_linesplit[0], straing_linesplit[0])
                rna_both_dict["different"].append(outstr)
        else:
            rna_ref_only_list.append(rna_dicts_list[0][ref_genename])
    
    for strain_genename in rna_strain_dicts_list[0]:
        if strain_genename not in rna_dicts_list[0]:
            rna_strain_only_list.append(rna_strain_dicts_list[0][strain_genename])

    #tRNA and rRNA
    for i in range(1,3):
        for ref_des in rna_dicts_list[i]:
            refgene_list = rna_dicts_list[i][ref_des]
            if ref_des in rna_strain_dicts_list[i]:
                straingene_list = rna_strain_dicts_list[i][ref_des]
                #easiest case
                if len(refgene_list) == 1 and len(straingene_list) == 1:
                    refg_linesplit = refgene_list[0].split("\t")
                    straing_linesplit = straingene_list[0].split("\t")
                    refg_len = refg_linesplit[2]
                    straing_len = straing_linesplit[2]
                    if refg_len == straing_len:
                        #description length id_ref id_strain
                        outstr = "{0}\t{1}\t{2}\t{3}\n".format(refg_linesplit[-1].strip(), refg_len, refg_linesplit[0], straing_linesplit[0])
                        rna_both_dict["same"].append(outstr)
                    else:
                        #description length_ref length_strain id_ref id_strain
                        outstr = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(refg_linesplit[-1].strip(), refg_len, straing_len, refg_linesplit[0], straing_linesplit[0])
                        rna_both_dict["different"].append(outstr)
                else:
                    for refgene in refgene_list:
                        uncommon_list.append(refgene)
                    for straingene in straingene_list:
                        uncommon_list.append(straingene)
            else:
                if len(refgene_list) > 1:
                    for element in refgene_list:
                        uncommon_list.append(element)
                    print("take a look! only in ref t/rRNA list longer 1")
                rna_ref_only_list.append(refgene_list[0])

        for strain_des in rna_strain_dicts_list[i]:
            if strain_des not in rna_dicts_list[i]:
                if len(rna_strain_dicts_list[i][strain_des]) > 1:
                    for element in refgene_list:
                        uncommon_list.append(element)
                    print("take a look! only in strain t/rRNA list longer 1")
                rna_strain_only_list.append(rna_strain_dicts_list[i][strain_des][0])
    
    print("#same RNAs: " + str(len(rna_both_dict["same"])))
    rna_both_same_path = f"{outpath_strain}_rna_both_same.tsv"
    header = "gene\tproduct\tref_length_bp\tref_locus_tag\tstrain_locus_tag\n"
    write_out_comparison(rna_both_same_path, rna_both_dict["same"], header)

    print("#RNAs only in strain: " + str(len(rna_strain_only_list)))
    print("#RNAs only in reference: " + str(len(rna_ref_only_list)))
    if len(rna_strain_only_list) >=1 or len(rna_ref_only_list) >=1:
        rna_only_in1_path = f"{outpath_strain}_rna_only_in1.tsv"
        rna_only_in1_file = open(rna_only_in1_path, "w")
        rna_only_in1_file.write("#only in strain")
        for entry in rna_strain_only_list:
            rna_only_in1_file.write(entry)
        rna_only_in1_file.write("\n#only in reference")
        for entry in rna_ref_only_list:
            rna_only_in1_file.write(entry)
        rna_only_in1_file.close()

    print("#different RNAs: " + str(len(rna_both_dict["different"])))
    print("#uncommon RNAs: " + str(len(uncommon_list)))
    rna_both_diff_path = f"{outpath_strain}_rna_both_diff.tsv"
    rna_both_diff_file = open(rna_both_diff_path, "w")
    for entry in rna_both_dict["different"]:
        rna_both_diff_file.write(entry)
    rna_both_diff_file.write("\n")
    for entry in uncommon_list:
        rna_both_diff_file.write(entry)
    rna_both_diff_file.close()


#prokka abgleich reference
# first gene abgleich, dann noch cog product abgleichen
# länge der gene mit aufnehmen
#output: genes in both list (tab separated) gene description len_in_ref len_in_strain
#output: genes only in one

ref_inpath = "results/analysis/ref_prokka/"
strain_inpath = "results/analysis/prokka/"
outpath = "results/analysis/ref_vs_prokka/"

    #"GCF_ecloacae", "GCF_abaumannii", "GCF_kpneumoniae", "GCF_ecolik12"
ref_name = "GCF_abaumannii"
os.system(f"mkdir {outpath}{ref_name}")
filter_hypothetical_proteins(f"{ref_inpath}{ref_name}/{ref_name}.tsv")

    #"13", "62", "66", "102", "109", "115", "129", "167", "169", "265", "281"
strain = "115"
os.system(f"mkdir {outpath}{strain}")
filter_hypothetical_proteins(f"{strain_inpath}{strain}/{strain}.tsv")

outpath_ref = f"{outpath}{ref_name}/{ref_name}"
outpath_strain = f"{outpath}{strain}/{strain}"

ref_tsv_path = f"{outpath_ref}_wo_hyprot.tsv"
seperate_ftypes(outpath, ref_tsv_path, ref_name)
rna_tsv_path = f"{outpath_ref}_rna.tsv"
cds_tsv_path = f"{outpath_ref}_cds.tsv"
cds_wog_tsv_path = f"{outpath_ref}_cds_wog.tsv"

    #[0]: tmRNA, [1]: tRNA, [2]: rRNA
rna_dicts_list = make_rna_dicts(rna_tsv_path)
cds_dict = make_cds_dict(cds_tsv_path)
    #[cds_wog_cogs_dict, cds_wog_ec_dict, cds_wog_only_des_dict]
cds_wog_dict_list = make_cds_wog_dicts(cds_wog_tsv_path)

strain_tsv_path = f"{outpath_strain}_wo_hyprot.tsv"
seperate_ftypes(outpath, strain_tsv_path, strain)
strain_rna_tsv_path = f"{outpath_strain}_rna.tsv"
strain_cds_tsv_path = f"{outpath_strain}_cds.tsv"
cds_wog_strain_tsv_path = f"{outpath_strain}_cds_wog.tsv"

    #[0]: tmRNA, [1]: tRNA, [2]: rRNA
rna_strain_dicts_list = make_rna_dicts(strain_rna_tsv_path)
rna_ref_vs_strain(rna_dicts_list, rna_strain_dicts_list, outpath_strain)

strain_cds_dict = make_cds_dict(strain_cds_tsv_path)
cds_ref_vs_strain(cds_dict, strain_cds_dict, outpath_strain)

    #[cds_wog_cogs_dict, cds_wog_ec_dict, cds_wog_only_des_dict]
cds_wog_strain_dict_list = make_cds_wog_dicts(cds_wog_strain_tsv_path)
cds_wog_ref_vs_strain(cds_wog_dict_list, cds_wog_strain_dict_list, outpath_strain)

"""
#old methods
def make_rna_ref_vs_strain_dicts(strain_tsv_path):
    strain_tsv_file = open(strain_tsv_path, "r")
    # dicts for reference strain comparision
    trna_rs_dict = {"both":[], "lonly":[]}
    rrna_rs_dict = {"both":[], "lonly":[]}
    tmrna_rs_dict = {"both":[], "lonly":[]}

    for gene in strain_tsv_file:
        splitted = gene.split("\t")
        ftype = splitted[1]
        gene_name = splitted[3]
        description = splitted[-1].strip()
        if ftype == "tmRNA":
            if gene_name in rna_dicts_list[0]:
                refgene = rna_dicts_list[0][gene_name]
                ref_id = refgene.split("\t")[0]
                ref_len = refgene.split("\t")[2]

                outstr = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(gene_name, description, splitted[2], ref_len, splitted[0], ref_id)
                tmrna_rs_dict["both"].append(outstr)
            else:
                outstr = "{0}\t{1}\t{2}\t{3}\n".format(gene_name, description, splitted[2], splitted[0])
                tmrna_rs_dict["lonly"].append(outstr)

        if ftype == "tRNA":
            if description in rna_dicts_list[1]:
                refgene_list = rna_dicts_list[1][description]
                if len(refgene_list) == 1:
                    ref_id = refgene_list[0].split("\t")[0]
                    ref_len = refgene_list[0].split("\t")[2]
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(description, splitted[2], ref_len, splitted[0], ref_id)
                else:
                    ref_outstr = ""
                    for refgene in refgene_list:
                        ref_id = refgene.split("\t")[0]
                        ref_len = refgene.split("\t")[2]
                        ref_outstr = ref_outstr + "{0}\t{1}\n".format(ref_len, ref_id)
                    outstr = "{0}\t{1}\t{2}\n{3}".format(description, splitted[2], splitted[0], ref_outstr)

                trna_rs_dict["both"].append(outstr)
            else:
                outstr = "{0}\t{1}\t{2}\n".format(description, splitted[2], splitted[0])
                trna_rs_dict["lonly"].append(outstr)
        
        if ftype == "rRNA":
            if description in rna_dicts_list[2]:
                refgene_list = rna_dicts_list[2][description]
                
                if len(refgene_list) == 1:
                    ref_id = refgene_list[0].split("\t")[0]
                    ref_len = refgene_list[0].split("\t")[2]
                    outstr = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(description, splitted[2], ref_len, splitted[0], ref_id)
                    if len(refgene_list) == 1:
                        outstr = "{0}\t{1}\t{2}\t{3}\t{4}\n".format(description, splitted[2], ref_len, splitted[0], ref_id)
                else:
                    ref_outstr = ""
                    for refgene in refgene_list:
                        ref_id = refgene.split("\t")[0]
                        ref_len = refgene.split("\t")[2]
                        ref_outstr = ref_outstr + "{0}\t{1}\n".format(ref_len, ref_id)
                    outstr = "{0}\t{1}\t{2}\n{3}".format(description, splitted[2], splitted[0], ref_outstr)

                rrna_rs_dict["both"].append(outstr)
            else:
                outstr = "{0}\t{1}\t{2}\n".format(description, splitted[2], splitted[0])
                rrna_rs_dict["lonly"].append(outstr)
    
    return([tmrna_rs_dict, trna_rs_dict, rrna_rs_dict])

def write_rna_outfiles(rna_rs_dicts_list, outpath_strain):
    frna_both_path = f"{outpath_strain}_rna_both.tsv"
    frna_both = open(frna_both_path, "w")
    frna_lonly_path = f"{outpath_strain}_rna_lonly.tsv"
    frna_lonly = open(frna_lonly_path, "w")

    for dict in rna_rs_dicts_list:
        both_list = dict["both"]
        print(len(both_list))
        if len(both_list) > 1:
            both_list.sort()
            for value in both_list:
                frna_both.write(value)
        else:
            frna_both.write(both_list[0])
        frna_both.write("\n")

        lonly_list = dict["lonly"]
        print(len(lonly_list))
        if len(lonly_list) == 0:
            continue
        elif len(lonly_list) >= 1:
            lonly_list.sort()
            for value in lonly_list:
                frna_lonly.write(value)
        elif len(lonly_list) == 1:
            frna_lonly.write(lonly_list[0])
        frna_lonly.write("\n")
    frna_both.close()
    frna_lonly.close()



#statt händisch:

frna_both = open("results/analysis/ref_vs_prokka/13/rna_both.tsv", "r")
new_out = []
des = ""
for line in frna_both:
    if line == "\n":
        new_out.append(line)
        continue
    splitted = line.split("\t")
    if len(splitted) > 3:
        is_double = False
        new_out.append(line)
    elif len(splitted) == 3:
        last_des = des
        is_double = True
        des = splitted[0]
    else:
        is_double = True

    if is_double:


   
        tmrna_dict[gene_name] = gene
    elif ftype == "tRNA":
        if description in trna_dict:
            trna_dict[description].append(gene)
        else:
            trna_dict[description] = [gene]
    elif ftype == "rRNA":
        if description in rrna_dict:
            rrna_dict[description].append(gene)
        else:
            rrna_dict[description] = [gene]

des_len_id = value.split("\n")
            des_len = des_len_id[:des_len_id.rfind("\t")]
            if des_len not in des_len_list:
                des_len_list.append(des_len)
                number = 1
            else:
                number += 1
"""
"""
trna_dict = {}
rrna_dict = {}
tmrna_dict = {}
cds_dict = {}
cds_wo_genename = []
genes_in_ref = []
for gene in tsv:
    splitted = gene.split("\t")
    ftype = splitted[1]
    gene_name = splitted[3]
    #len [0] and description [1] (#gene_len = splitted[2])
    description = splitted[-1].strip()
    len_gene_list = [splitted[2], gene]

    if ftype == "CDS":
        number_cds += 1
        if gene_name == "":
            cds_wo_genename.append(gene)
            genes_in_ref.append(description)
        elif gene_name in cds_dict:
            cds_dict[gene_name].append(len_gene_list)
        else:
            genes_in_ref.append(gene_name)
            cds_dict[gene_name] = [len_gene_list]

    elif ftype == "tmRNA":
        number_tmrna += 1
        if gene_name in tmrna_dict:
            tmrna_dict[gene_name].append(len_gene_list)
        else:
            tmrna_dict[gene_name] = [len_gene_list]

    elif ftype == "tRNA":
        number_trna += 1
        if description in trna_dict:
            trna_dict[description].append(len)
        else:
            trna_dict[description] = [len]

    elif ftype == "rRNA":
        number_rrna += 1
        if description in rrna_dict:
            rrna_dict[description].append(len)
        else:
            rrna_dict[description] = [len]
tsv.close()

print("CDS without genename: " + str(len(cds_wo_genename)))
print("CDS in reference: " + str(len(genes_in_ref)))
print(f"#CDS: {number_cds}")
print(f"#tmRNA: {number_tmrna}")
print(f"#tRNA: {number_trna}")
print(f"#rRNA: {number_rrna}")

#genes_in_both = "results/analysis/ref_vs_prokka/13/both.tsv"
#genes_in_strain = "results/analysis/ref_vs_prokka/13/only_strain.tsv"
#genes_in_ref = "results/analysis/ref_vs_prokka/13/only_ref.tsv"

cds_both_tsv = open("results/analysis/ref_vs_prokka/13/cds_both.tsv", "w")
cds_only_strain_tsv = open("results/analysis/ref_vs_prokka/13/cds_only_strain.tsv", "w")
cds_only_ref_tsv = open("results/analysis/ref_vs_prokka/13/cds_only_ref.tsv", "w")

strain_path = "results/analysis/prokka/13/13_wo_hyprot.tsv"
fstrain = open(strain_path, "r")
genes_in_strain = []

number_only_strain = number_only_ref = number_both = 0
for gene in fstrain:
    splitted = gene.split("\t")
    ftype = splitted[1]
    gene_name = splitted[3]
    gene_id = splitted[0]
    #len [0] and description [1] (#gene_len = splitted[2])
    gene_len = splitted[2]
    description = splitted[-1].strip()
    in_both = False

    if ftype == "CDS":
        if gene_name in cds_dict:
            in_both = True
            genes_in_ref.remove(gene_name)
            if len(cds_dict[gene_name]) > 1:
                print("yes longer 1")
                ref_len = []
                ref_id = []
                for refgene in cds_dict[gene_name]:
                    id = refgene[1].split("\t")[0]
                    ref_id.append(id)
                    ref_len.append(refgene[0])
            else:
                ref_len = cds_dict[gene_name][0][0]
                ref_line = cds_dict[gene_name][0][1]
                ref_id = ref_line.split("\t")[0]
            both = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(gene_name, description, gene_len, ref_len, gene_id, ref_id)
        """"""
        else:
            for uncom_gene in cds_wo_genename:
                if uncom_gene.find(description) >=0:
                    in_both = True
                    if description in genes_in_ref:
                        genes_in_ref.remove(description)
                    ref_len = uncom_gene.split("\t")[2]
                    both = "{0}\t{1}\t{2}\t{3}\n".format(gene_name, description, gene_len, ref_len)
                    break
                    """"""
        if in_both:
            number_both +=1
            cds_both_tsv.write(both)
        else:
            number_only_strain +=1
            cds_only_strain_tsv.write(gene)

if len(genes_in_ref) > 0:
    for lonly_gene in genes_in_ref:
        if lonly_gene in cds_dict:
            number_only_ref +=1
            cds_only_ref_tsv.write(cds_dict[lonly_gene][0][1])
        """"""
        else:
            for uncom_gene in cds_wo_genename:
                if uncom_gene.find(lonly_gene) >=0:
                    number_only_ref +=1
                    cds_only_ref_tsv.write(uncom_gene)
                    break
                """"""

cds_both_tsv.close()
cds_only_strain_tsv.close()
cds_only_ref_tsv.close()

print(f"#genes in both: {number_both}")
print(f"#genes only in strain: {number_only_strain}")
print(f"#genes only in reference: {number_only_ref}")
"""