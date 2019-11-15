"""
for i in `ls ./*.pbs`; do qsub  $i;done
"""

import os
import pandas as pd
import re

def create_pbs():
    files = os.listdir("/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/pooled_peak_files")
    
    for one_file in files:
        name = one_file.split("_")[0]
        handle = open(name + ".pbs", "w")
        file_path =os.path.join("/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/pooled_peak_files", one_file)
        output_dir = "/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/output/%s" % name
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        handle.write("findMotifsGenome.pl %s mm10 %s" %(file_path, output_dir))
        handle.close()
        
def annot_promoter():
    files = os.listdir("/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/pooled_peak_files")
    
    for one_file in files:
        name = one_file.split("_")[0]
        handle = open(name + ".annotatePeaks.pbs", "w")
        file_path =os.path.join("/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/pooled_peak_files", one_file)
        output_file = "/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/pooled_peak_files_remove_promoter/%s.bed" % name
        handle.write("annotatePeaks.pl  %s mm10 > %s" %(file_path, output_file))
        handle.close()
    
def remove_promoter():
    annotation_dir = "/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/pooled_peak_files_remove_promoter"
    columns = ['PeakID',
       'Chr', 'Start', 'End', 'Strand', 'Peak Score',
       'Focus Ratio/Region Size', 'Annotation', 'Detailed Annotation',
       'Distance to TSS', 'Nearest PromoterID', 'Entrez ID', 'Nearest Unigene',
       'Nearest Refseq', 'Nearest Ensembl', 'Gene Name', 'Gene Alias',
       'Gene Description', 'Gene Type']
    #df_tmp = pd.DataFrame(columns=columns)
    result_file_list = os.listdir(annotation_dir)
    
    for one_file in result_file_list:
        if re.search("rmpro", one_file):
            continue
        file_path = os.path.join(annotation_dir, one_file)
        df = pd.read_csv(file_path, sep="\t", header=0)
        df.columns = columns
        df_tmp = df[(df["Distance to TSS"] > 1500) | (df["Distance to TSS"] < -1500) ]
        df_tmp = df_tmp[['PeakID', 'Chr', 'Start', 'End', 'Strand', 'Peak Score']]
        df_tmp.to_csv(file_path + ".rmpro", header=False, index=False, sep="\t")
        #exit(0)
        
def create_new_pbs():
    files = os.listdir("/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/pooled_peak_files_remove_promoter")
    
    for one_file in files:
        if not re.search("rmpro", one_file):
            continue
        name = one_file.split(".")[0]
        handle = open(name + ".pbs", "w")
        file_path =os.path.join("/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/pooled_peak_files_remove_promoter", one_file)
        output_dir = "/home/zhluo/Project/CRC/data_nazhang/step34_TF_enrichment/output_new/%s" % name
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        handle.write("findMotifsGenome.pl %s mm10 %s" %(file_path, output_dir))
        handle.close()
       
 
if __name__ == "__main__":    
    #annot_promoter()
    remove_promoter()
    create_new_pbs()
