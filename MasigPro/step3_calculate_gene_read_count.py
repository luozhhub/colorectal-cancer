import sys,re, os
import pandas as pd
import numpy as np

class cal_gene_marks():
    def __init__(self):
        self.accession_path = "/home/zhluo/Project/CRC/TF_enrichment/homer/data/accession/mouse2gene.tsv"
        
        
    def parse_ID(self):
        df = pd.read_csv(self.accession_path, names=["complex", "entrez", "unigene", "RefSeq", "Ensembl", "symbol"], sep="\t", index_col=False)
        #print(df[0:5])
        df_sub = df[["RefSeq", "Ensembl"]]
        #df_sub =df_sub[(df_sub["RefSeq"] != np.nan) & (df_sub["Ensembl"] != np.nan)]
        df_sub = df_sub.dropna()
        df_sub = df_sub.drop_duplicates()
        print(len(df_sub))
        
    def get_gene_body(self):
        df = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/refernece_genome/genecode/gene_locus.bed", names=["chr", "start", "end", "strand", "ensembl"], sep="\t", index_col=False)
        ensembl = [item.strip("\"").split(".")[0] for item in df["ensembl"]]
        df["ensembl"] = ensembl
        print(len(df))
        print(len(set(list(df["ensembl"]))))
        df.to_csv("/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_body.bed", sep="\t", header=False,  index=False)
        

if __name__ == "__main__":
    gene_cal = cal_gene_marks()
    #gene_cal.parse_ID()
    gene_cal.get_gene_body()