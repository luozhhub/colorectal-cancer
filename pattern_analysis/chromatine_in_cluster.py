import sys,re, os
import pandas as pd
import numpy as np
import collections
import tempfile
import subprocess
from io import StringIO
from functools import reduce



class makeStateMatrix():
    def __init__(self):
        self.accession_path = "/home/zhluo/Project/CRC/TF_enrichment/homer/data/accession/mouse2gene.tsv"
        
    def writeToTmp(self, dataframe):
        tempFd, tempPath = tempfile.mkstemp()
        dataframe.to_csv(tempPath, sep="\t", header=False, index=False)
        os.close(tempFd)
        return tempPath
        
    def get_gene_body(self, extend=10000):
        df = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/refernece_genome/genecode/gene_locus.bed", \
                          names=["chr", "start", "end", "strand", "ensembl"], sep="\t", index_col=False)
        ensembl = [item.strip("\"").split(".")[0] for item in df["ensembl"]]
        df["ensembl"] = ensembl
        df["start"] = [item - extend for item in df["start"]]
        df["end"] = [item + extend for item in df["end"]]
        df.loc[df["start"] < 0, "start"] = 0
        print(len(df))
        print(len(set(list(df["ensembl"]))))
        bed_file = self.writeToTmp(df)
        return bed_file
        
    def intersect_with_state(self, bed_file):
        chromstate_dir = "/home/zyang/Project/CRC/step45_chromatin_state_in_cluster/state13"
        df_list = []
        columns = []
        for one_file in os.listdir(chromstate_dir):
            
            if "segments.bed" in one_file:
                df = pd.read_csv(os.path.join(chromstate_dir, one_file), sep="\t", names=["chr", "start", "end", "state"])
                df = df[df.state != "E2"]
                df_tmp = self.writeToTmp(df)
                #tmp_result = "/home/zyang/Project/CRC/step45_chromatin_state_in_cluster/test.txt"
                cmd = "bedtools intersect -wao -a %s -b %s" % (bed_file, df_tmp)
                out_string = subprocess.check_output([cmd], shell=True)
                out_string = StringIO(out_string.decode('utf-8'))
                df_result = pd.read_csv(out_string, sep="\t", \
                            names=["chr", "start", "end", "strand", "ensembl", "chr_state", "start_state", "end_state", "state", "overlap"])
                df_result = df_result[["ensembl", "state", "overlap"]]
                
                df_result = df_result[df_result.state != "."].reset_index()
                df_result = df_result.groupby(['ensembl','state'])["overlap"].sum().reset_index()
                
                idx = df_result.groupby(['ensembl'])['overlap'].transform(max) == df_result['overlap']
                df_result = df_result[idx]
                counts = collections.Counter(df_result.ensembl)
                ##some gene has several equit max value
                del_ensembl = []
                for key in counts:
                    if counts[key] > 1:
                        del_ensembl.append(key)
                df_result = df_result[~df_result["ensembl"].isin(del_ensembl)]        
                print(len(set(df_result.ensembl)))
                print(len(df_result.ensembl))
                #df_result.to_csv("test1.txt", sep="\t", index=True)
                #os.remove(tmp_result)
                os.remove(df_tmp)
                df_result = df_result[["ensembl", "state"]]
                df_list.append(df_result)
                columns.append(one_file.split("_")[0])
                #print(df_result[0:5])
                
            
        df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['ensembl'], how='outer'), df_list)
        df_merged.index = df_merged["ensembl"]
        df_merged.drop(['ensembl'], axis=1, inplace=True)
        df_merged.columns = columns
        df_merged = df_merged[["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]]
        #df_merged.to_csv("state_summary_cluster.txt", sep="\t", header=True)
        return df_merged
                
                  #return 0
    def gene_list_in_RNA_cluster(self):
        RNA_gene_list_dir = "/home/zyang/Project/CRC/step40_masigpro_gene_list/rna_gene_list"
        dict_geneList = {}
        for one_file in os.listdir(RNA_gene_list_dir):
            if "rna_gene_list_cluster" in one_file:
                df_gene_cluster = pd.read_csv(os.path.join(RNA_gene_list_dir, one_file), sep= "\t")
                gene_list = list(df_gene_cluster["ensembl"])
                name = one_file.split(".")[0]
                dict_geneList[name] = gene_list
        return dict_geneList
        
    def main(self):
        bed_file = self.get_gene_body()
        print(bed_file)
        df_merged = self.intersect_with_state(bed_file)
        #df_merged = pd.read_csv("/home/zyang/Project/CRC/step45_chromatin_state_in_cluster/state_summary_cluster.txt", sep="\t", index_col=0)
        os.remove(bed_file)
        dict_geneList = self.gene_list_in_RNA_cluster()
        sub_df = pd.DataFrame(columns=["ctrl", "2weeks", "4weeks", "7weeks", "10weeks", "cluster"])
        
        for key in dict_geneList:
            gene_list = dict_geneList[key]
            #print(gene_list)
            #return 0
            tmp_df = df_merged[df_merged.index.isin(gene_list)]
            tmp_df["cluster"] = key
            sub_df = sub_df.append(tmp_df)
        sub_df.to_csv("state_summary_cluster.txt", sep="\t", header=True)
            
            
        
if __name__ == "__main__":
    mkmatrix = makeStateMatrix()
    mkmatrix.main()
    
    