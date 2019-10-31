"""
# in workstation :/home/zhihl/Project/ZNF_database/cistrome/mouse_bed/mouse_factor/
for i in `awk -F"\t" '{print $8}' ./RELA_df.txt`; do scp "/home/zhihl/Project/ZNF_database/cistrome/mouse_bed/mouse_factor/"$i zhluo@xxxx:/home/zhluo/Project/CRC/data_nazhang/step36_NFkB/RELA_data/ ;done
#eg.
annotatePeaks.pl ./RELA_data/5183_sort_peaks.narrowPeak.bed mm10 > ./annotatePeaks_output/5183_sort_peaks.narrowPeak.bed.anno.txt
#run in cluster
cd ./RELA_data
for i in `ls ./` ;do annotatePeaks.pl $i mm10 > "../annotatePeaks_output/"${i%.*}".txt" ;done


after run this script
cd ~/Project/CRC/Chip_analysis/dff/version0821/final_result/target_gene
scp zhluo@xxxx:/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/TF_enrichment/gene_numebr.txt ./
after run select peak in step11_fpkm.Rmd
scp ./whole_result.bed zhluo@xxxx:/home/zhluo/Project/CRC/data_nazhang/step36_NFkB
qsub ./whole_H3K27ac.pbs
"""

import os,sys
import pandas as pd
from collections import Counter

def merge_df():
    annotation_dir = "/home/zhluo/Project/CRC/data_nazhang/step36_NFkB/annotatePeaks_output"
    columns = ['PeakID',
       'Chr', 'Start', 'End', 'Strand', 'Peak Score',
       'Focus Ratio/Region Size', 'Annotation', 'Detailed Annotation',
       'Distance to TSS', 'Nearest PromoterID', 'Entrez ID', 'Nearest Unigene',
       'Nearest Refseq', 'Nearest Ensembl', 'Gene Name', 'Gene Alias',
       'Gene Description', 'Gene Type']
    df_tmp = pd.DataFrame(columns=columns)
    result_file_list = os.listdir(annotation_dir)
    
    for one_file in result_file_list:
        file_path = os.path.join(annotation_dir, one_file)
        df = pd.read_csv(file_path, sep="\t", header=0)
        df.columns = columns
        df_tmp = df_tmp.append(df)
        
    gene_name_list = list(df_tmp["Nearest Ensembl"])
    gene_number = Counter(gene_name_list)
    return(gene_number)   
    
     
        
    
if __name__ == "__main__":
    gene_number_dict = merge_df()
    df_gene_number = pd.DataFrame(columns=["gene_name", "num"])
    n = 0
    for i in gene_number_dict:
        df_gene_number.loc[str(n)] = [i, gene_number_dict[i]]
        n = n + 1
    print(df_gene_number)
    df_gene = df_gene_number.sort_values(by=["num"], axis=0, na_position='last', ascending=False)
    df_gene = df_gene.reset_index(drop=True)
    df_gene.to_csv("./gene_numebr.txt", header=True, index=False, sep="\t")
    
    