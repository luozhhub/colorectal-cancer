import sys,re, os
import pandas as pd
import numpy as np
import collections

class cal_gene_marks():
    def __init__(self):
        self.accession_path = "/home/zhluo/Project/CRC/TF_enrichment/homer/data/accession/mouse2gene.tsv"
        
    def run_multiBamSummary(self, markers=["H3K27ac", "H3K4me1", "Input"], bs="enhancer_neargene.bed", bed_file="/home/zhluo/Project/CRC/data_nazhang/step32_element/enhancer/enhancer.peak.unique.ID.bed"):
        #calculate read count
        #bam files in /home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles
        bamDir = "/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles"
        output_dir = "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/multiBamSummary_gene_body"
        #bed_file = "/home/zhluo/Project/CRC/data_nazhang/step32_element/enhancer/enhancer.peak.unique.ID.bed"
        weeks = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
        replicate = ["1", "2" ,"3"]
        markers = markers
        suffix = "_combined_1_paired.merged.nodup.bam"
        sample_list = []    
        lable_list = []
        for we in weeks:
            for mark in markers:
                for rep in replicate:
                    sample = "%s-%s-%s%s" %(we, rep, mark, suffix)
                    sample_list.append(os.path.join(bamDir, sample))
                    lable_list.append("%s_%s_%s" % (we, rep, mark)) 
        bamfiles = " ".join(sample_list)
        lables = " ".join(lable_list)
        cmd = "multiBamSummary BED-file --BED %s --numberOfProcessors 25 --bamfiles %s  --minMappingQuality 30 --labels %s -out %s --outRawCounts %s" %(bed_file, bamfiles, \
        lables, os.path.join(output_dir, "readCounts_%s.npz" % bs), os.path.join(output_dir, "readCounts_%s.tab" % bs))
        pbs_handle = open("/home/zhluo/Project/CRC/data_nazhang/step39_read_count/run_multiBamSummary_%s.pbs" % bs, "w")
        pbs_handle.write(cmd)
        pbs_handle.close()
        
        
    def merge_peakID(self):
        
        marker = "enhancer_neargene"
        enhancer = "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/multiBamSummary_gene_body/readCounts_enhancer_neargene.bed.tab"
        enhancer_gene = "/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/enhancer.peak.unique.ID.gene_name.bed"
        read_count_table = enhancer
        
        gene_table = enhancer_gene
        #chr     start   end     ID      middle  gene_id
        gene_df = pd.read_csv(gene_table, sep="\t", index_col=False)
        gene_df.columns = ["chr", "start" ,"end",  "peakID", "middle",  "ensembl"]
        #gene_df.index = gene_df["peakID"]
        
        gene_df = gene_df[["chr", "start", "end", "ensembl"]]
        print(gene_df[0:5])
        #exit(1)
        marker_df = pd.read_csv(read_count_table, sep="\t", index_col=False, header=0, quoting=3)
        marker_df.columns = [item.strip("#'") for item in marker_df]
        marker_merge_df = pd.merge(marker_df, gene_df,  how='left', left_on=["chr", "start", "end"], right_on = ["chr", "start", "end"])
        print(marker_merge_df[0:5])
        marker_merge_df = marker_merge_df.drop_duplicates()
        marker_merge_df = marker_merge_df.dropna()
        marker_merge_df.index = marker_merge_df["ensembl"]
        marker_merge_df = marker_merge_df.drop(["chr", "start", "end", "ensembl"], axis=1)
        marker_merge_df = marker_merge_df.groupby(marker_merge_df.index).sum()
                
        #pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/enhancer.peak.unique.ID.gene_name.bed", sep="\t", index_col=False, names=["chr", "start", "end", "peakID"])
        marker_merge_df.to_csv("/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/new_bed_for_bw/enhancer_total_readCount.txt" , sep="\t", header=True,  index=True)
        
        
    
        
if __name__ == "__main__":
    cal_gene = cal_gene_marks()
    #cal_gene.run_multiBamSummary()
    """
    qsub -q batch -V -l nodes=1:ppn=26 ./run_multiBamSummary_enhancer_neargene.bed.pbs
    """
    cal_gene.merge_peakID()