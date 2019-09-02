"""
this script is used to find snp nearest gene
"""
#!/usr/bin/pyhton
import pandas as pd
import json

def read_LD():
    return json.load(open('./ld_dict.json', 'r'))

def find_coding_gene():
    """
    rsID ID_nearest_gene_snpsnap ID_nearest_gene_snpsnap_protein_coding HGNC_nearest_gene_snpsnap
    """
    annotation_file = "/home/zhluo/Project/CRC/LD_data/ensembl_io/snpsnap_gene.txt"
    df = pd.read_csv(annotation_file, sep="\t", header=0)
    print (df.head())
    id_dict = {}
    
    for index, row in df.iterrows():
        id = row["rsID"]
        if not id.startswith("rs"):
            continue
        gene = row["ID_nearest_gene_snpsnap"]
        id_dict[id] = gene
    
     
    snp_dict = read_LD()
    whole_list = []
    for key in snp_dict.keys():
        snp_list = snp_dict[key] + [key]
        whole_list = whole_list + snp_list
    whole_list = list(set(whole_list))
    output = open("snp_gene_map.txt", "w")
    for one_rs in whole_list:
        if one_rs in id_dict:
            output.write("%s\t%s\n" %(one_rs, id_dict[one_rs]))
        else:
            output.write("%s\t%s\n" %(one_rs, None))
    output.close()
    #or line in open(annotation_file, "r"):
    
        
if __name__ == "__main__":
    find_coding_gene()

