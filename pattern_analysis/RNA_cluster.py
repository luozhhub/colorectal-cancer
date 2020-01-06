import os,sys
import pandas as pd

def DEG_in_cluster():
    dict_geneList = {}
    RNA_gene_list_dir = "/home/zyang/Project/CRC/step40_masigpro_gene_list/rna_gene_list"
    step_dir = "/home/zyang/Project/CRC/step42_TF_enrichment/"
    
    for one_file in os.listdir(RNA_gene_list_dir):
        if "rna_gene_list_cluster" in one_file:
            df_gene_cluster = pd.read_csv(os.path.join(RNA_gene_list_dir, one_file), sep= "\t")
            gene_list = list(df_gene_cluster["ensembl"])
            name = one_file.split(".")[0]
            dict_geneList[name] = gene_list
            
    df_total = pd.DataFrame(columns=["cluster", "marker", "week", "deg_number"])
    DEG_gene_dir = "/home/zyang/Project/CRC/step43_fold_change_gene"
    for one_file in os.listdir(DEG_gene_dir):
        if ".bed" in one_file:
            [marker, week] = one_file.strip(".bed").split("_")[1:3]
            #print(marker, week)
            df_deg = pd.read_csv(os.path.join(DEG_gene_dir, one_file), sep="\t", names=["chr", "start", "end", "peak_id", "middle", "ensembl"])
            deg_genes = list(set(list(df_deg["ensembl"])))
            for key in dict_geneList:
                gene_list = dict_geneList[key]
                number = len(list(set(deg_genes) & set(gene_list)))
                df_total = df_total.append(pd.Series({"cluster": key, "marker": marker, "week": week.strip("weeks"), "deg_number": number}), ignore_index=True)
            
    output_dir = "/home/zyang/Project/CRC/step43_fold_change_gene/rna_cluster"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    df_total.to_csv(os.path.join(output_dir, "DEG_summary.txt"), sep="\t", header=True)
            
        
        
        
        
        
        
        
        

if __name__ == "__main__":
    DEG_in_cluster()