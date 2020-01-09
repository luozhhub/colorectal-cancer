import os,sys
import pandas as pd

def translate_human_symbol_ensembl(symbol_list):
    mart_table = pd.read_csv("/home/zyang/Project/CRC/step44_driver_gene/mart_export_human_geneid_symbol.txt", sep="\t", quotechar='"')
    result_df = mart_table[mart_table["external_gene_name"].isin(symbol_list)]
    return list(set(list(result_df["ensembl_gene_id"])))
    
def translate_human_to_mouse(ensembl_list):
    mart_table = pd.read_csv("/home/zyang/Project/CRC/step44_driver_gene/mart_export_humanGene_mouseGene.txt", sep="\t", quotechar='"')
    mart_table.columns = ["human_gene", "human_transcript", "mouse_gene"]
    result_df = mart_table[mart_table["human_gene"].isin(ensembl_list)]
    return list(set(list(result_df["mouse_gene"])))

def gene_list_in_RNA_cluster():
    RNA_gene_list_dir = "/home/zyang/Project/CRC/step40_masigpro_gene_list/rna_gene_list"
    dict_geneList = {}
    for one_file in os.listdir(RNA_gene_list_dir):
        if "rna_gene_list_cluster" in one_file:
            df_gene_cluster = pd.read_csv(os.path.join(RNA_gene_list_dir, one_file), sep= "\t")
            gene_list = list(df_gene_cluster["ensembl"])
            name = one_file.split(".")[0]
            dict_geneList[name] = gene_list
    return dict_geneList

def parse_driver_gene():
    """
    data are download from driverDBv3: http://driverdb.tms.cmu.edu.tw/download
    """
    tab_file = open("/home/zyang/Project/CRC/step44_driver_gene/mutation_download_tab.txt", "r")
    all_data = tab_file.readlines()
    coad_tmp = []
    for line in all_data:
        array = line.strip().split("\t")
        cancer_type, method, gene_list = array[1], array[2], array[3].split(",")
        gene_list = [item.strip(" ") for item in gene_list]
        if cancer_type == "COAD":
            coad_tmp = coad_tmp + gene_list
    coad = list(set(coad_tmp))
    return coad
    
    
def main():
    colon_driver_gene = translate_human_to_mouse(translate_human_symbol_ensembl(parse_driver_gene()))
    gene_dict = gene_list_in_RNA_cluster()
    df_driver = pd.DataFrame(columns=["driver_gene_number", "cluster_name", "cluster_gene_number", "intersect_number"])
    driver_number = len(colon_driver_gene)
    for key in gene_dict:
        gene_list = gene_dict[key]
        number = len(list(set(colon_driver_gene) & set(gene_list)))
        
        df_driver = df_driver.append(pd.Series({"driver_gene_number" : driver_number, "cluster_name": key, \
        "cluster_gene_number" : len(gene_list), "intersect_number": number}), ignore_index=True)
    df_driver.to_csv("/home/zyang/Project/CRC/step44_driver_gene/driver_summary.txt", sep="\t", header=True)
        
if __name__ == "__main__":
    main()

    