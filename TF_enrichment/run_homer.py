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
        
#201200103 select RNA gene list to do TF enrichment
#the TSS file path: /home/zyang/Project/CRC/step42_TF_enrichment/mm10.Tss.gff
#for H3K27ac and H3K4me1
def select_peaks(step=1):
    all_peak_dir = "/home/zyang/Project/CRC/step42_TF_enrichment/all_peaks"
    if step <1:
        
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
            df_tmp = df_tmp[['PeakID', 'Chr', 'Start', 'End', 'Strand', 'Peak Score', 'Nearest Ensembl']]
            df_tmp.to_csv(os.path.join(all_peak_dir, one_file + ".ensembl.rmpro.txt"), header=False, index=False, sep="\t")
            
    if step <2:
        #select gene_list
        dict_geneList = {}
        RNA_gene_list_dir = "/home/zyang/Project/CRC/step40_masigpro_gene_list/rna_gene_list"
        step_dir = "/home/zyang/Project/CRC/step42_TF_enrichment/"
        
        for one_file in os.listdir(RNA_gene_list_dir):
            if "rna_gene_list_cluster" in one_file:
                df_gene_cluster = pd.read_csv(os.path.join(RNA_gene_list_dir, one_file), sep= "\t")
                gene_list = df_gene_cluster["ensembl"]
                name = one_file.split(".")[0]
                dict_geneList[name] = gene_list
        
        for one_file in os.listdir(all_peak_dir):
            for marker in ["H3K27ac", "H3K4me1"]:
                if marker in one_file:
                    prefix = one_file.split(".")[0]
                    subpeak_dir = os.path.join(step_dir, "rna", "rna_" + prefix)
                    print(subpeak_dir)
                    if not os.path.exists(subpeak_dir):
                        os.makedirs(subpeak_dir)
                    df_week = pd.read_csv(os.path.join(all_peak_dir, one_file), sep="\t", names=['PeakID', 'Chr', 'Start', 'End', 'Strand', 'Peak Score', 'Nearest Ensembl'])
                    for key in dict_geneList:
                        gene_list = dict_geneList[key]
                        df_cluster = df_week[df_week['Nearest Ensembl'].isin(gene_list)]
                        df_cluster.to_csv(os.path.join(subpeak_dir, key + ".ensembl.rmpro.txt"), header=False, index=False, sep="\t")
                        
def TF_enrichment():
    rna_peak_dir = "/home/zyang/Project/CRC/step42_TF_enrichment/rna/"
    for one_dir in os.listdir(rna_peak_dir):
        for one_file in os.listdir(os.path.join(rna_peak_dir, one_dir)):
            file_path = os.path.join(rna_peak_dir, one_dir, one_file)
            output_dir = "/home/zyang/Project/CRC/step42_TF_enrichment/rna_output/%s/%s" % (one_dir, one_file.split(".")[0])
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            handle = open("/home/zyang/Project/CRC/pbs/" + one_dir + one_file.split(".")[0] + ".pbs", "w")
            handle.write("/home/zhluo/Project/CRC/TF_enrichment/homer/bin/findMotifsGenome.pl %s mm10 %s" %(file_path, output_dir))
            handle.close()
            
                
                
                
                
       
 
if __name__ == "__main__":    
    #annot_promoter()
    #remove_promoter()
    #create_new_pbs()
    #select_peaks()
    TF_enrichment()
