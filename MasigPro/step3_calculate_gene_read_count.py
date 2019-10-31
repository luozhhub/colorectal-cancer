import sys,re, os
import pandas as pd
import numpy as np
import collections

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
        df["start"] = [item - 2000 for item in df["start"]]
        df["end"] = [item + 2000 for item in df["end"]]
        print(len(df))
        print(len(set(list(df["ensembl"]))))
        
        df.to_csv("/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_body_2k.bed", sep="\t", header=False,  index=False)
        
        
    def run_multiBamSummary(self, markers=["H3K27me3"], bs="H3K27me3.bed", bed_file=None):
        #calculate read count
        #bam files in /home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles
        bamDir = "/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles"
        output_dir = "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/multiBamSummary_gene_body"
        bed_file = bed_file
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
        
    def get_TSS(self, distence=2000):
        tss_file = "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/multiBamSummary_TSS_2kb/mm10.Tss.gff"
        df = pd.read_csv(tss_file, sep="\t", index_col=False, names=["chr", "tss", "strand", "gene_1", "gene_id", "gene_2", "transcript_id"])
        df_sub = df[["chr", "tss", "gene_id"]]
        df_sub["gene_id"] = [item.strip(";").split(".")[0] for  item in df_sub["gene_id"]]
        df_sub = df_sub.groupby("gene_id").agg({"chr": "first", "tss": [np.min , np.max ]})
        df_sub.columns = df_sub.columns.droplevel(0)
        df_sub.columns = ["chr", "start", "end" ]
        df_sub["ensembl"] = df_sub.index
        #print(df_sub.groupby("gene_id")["chr"].first())
        
        df_sub["start"] = [item - distence for item in df_sub["start"]]
        df_sub["end"] = [item + distence for item in df_sub["end"]]
        df_sub.to_csv("/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_tss_%s.bed" % distence, sep="\t", header=False,  index=False)
        
    def merge_read_count(self):
        enhancer = "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/multiBamSummary_gene_body/readCounts_enhancer.bed.tab"
        promoter = "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/multiBamSummary_gene_body/readCounts_H3K4me3.bed.tab" 
        reppressed = "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/multiBamSummary_gene_body/readCounts_H3K27me3.bed.tab"
        heterochromatin = "/home/zhluo/Project/CRC/data_nazhang/step39_read_count/multiBamSummary_gene_body/readCounts_K92_K93.bed.tab"
        
        enhancer_gene = "/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_tss_100000.bed"
        promoter_gene = "/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_tss_2000.bed"
        repressed_gene = "/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_body_2k.bed"
        heterochromatin_gene = "/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_tss_10000.bed"
        
        promoter_df = pd.read_csv(promoter, sep="\t", index_col=False, header=0, quoting=3)
        promoter_df.columns = [item.strip("#'") for item in promoter_df]
        promoter_gene_df = pd.read_csv(promoter_gene, sep="\t", index_col=False, names=["chr", "start", "end", "ensembl"])
        promoter_merge_df = pd.merge(promoter_df, promoter_gene_df,  how='left', left_on=["chr", "start", "end"], right_on = ["chr", "start", "end"])
        promoter_merge_df = promoter_merge_df.drop_duplicates()
        promoter_merge_df = promoter_merge_df.dropna()
        promoter_merge_df.index = promoter_merge_df["ensembl"]
        promoter_merge_df = promoter_merge_df.drop(["chr", "start", "end", "ensembl"], axis=1)
        
        enhancer_df = pd.read_csv(enhancer, sep="\t", index_col=False, header=0, quoting=3)
        enhancer_df.columns = [item.strip("#'") for item in enhancer_df]
        enhancer_gene_df = pd.read_csv(enhancer_gene, sep="\t", index_col=False, names=["chr", "start", "end", "ensembl"])
        enhancer_merge_df = pd.merge(enhancer_df, enhancer_gene_df,  how='left', left_on=["chr", "start", "end"], right_on = ["chr", "start", "end"])
        enhancer_merge_df = enhancer_merge_df.drop_duplicates()
        enhancer_merge_df = enhancer_merge_df.dropna()
        enhancer_merge_df.index = enhancer_merge_df["ensembl"]
        enhancer_merge_df = enhancer_merge_df.drop(["chr", "start", "end", "ensembl"], axis=1)
        
        total_df = pd.merge(enhancer_merge_df, promoter_merge_df,  how='outer', left_index=True, right_index=True)
        
        
        
        reppressed_df = pd.read_csv(reppressed, sep="\t", index_col=False, header=0, quoting=3)
        reppressed_df.columns = [item.strip("#'") for item in reppressed_df]
        reppressed_gene_df = pd.read_csv(repressed_gene, sep="\t", index_col=False, names=["chr", "start", "end", "strand", "ensembl"])
        reppressed_gene_df = reppressed_gene_df[["chr", "start", "end", "ensembl"]]
        reppressed_merge_df = pd.merge(reppressed_df, reppressed_gene_df,  how='left', left_on=["chr", "start", "end"], right_on = ["chr", "start", "end"])
        reppressed_merge_df = reppressed_merge_df.drop_duplicates()
        reppressed_merge_df = reppressed_merge_df.dropna()
        reppressed_merge_df.index = reppressed_merge_df["ensembl"]
        reppressed_merge_df = reppressed_merge_df.drop(["chr", "start", "end", "ensembl"], axis=1)
        
        total_df = pd.merge(total_df, reppressed_merge_df,  how='outer', left_index=True, right_index=True)
        #print(total_df.columns)
        #total_df = total_df.drop(["ensembl"], axis=1)
        
        
        
        heterochromatin_df = pd.read_csv(heterochromatin, sep="\t", index_col=False, header=0, quoting=3)
        heterochromatin_df.columns = [item.strip("#'") for item in heterochromatin_df]
        heterochromatin_gene_df = pd.read_csv(heterochromatin_gene, sep="\t", index_col=False, names=["chr", "start", "end", "ensembl"])
        heterochromatin_merge_df = pd.merge(heterochromatin_df, heterochromatin_gene_df,  how='left', left_on=["chr", "start", "end"], right_on = ["chr", "start", "end"])
        heterochromatin_merge_df = heterochromatin_merge_df.drop_duplicates()
        heterochromatin_merge_df = heterochromatin_merge_df.dropna()
        heterochromatin_merge_df.index = heterochromatin_merge_df["ensembl"]
        heterochromatin_merge_df = heterochromatin_merge_df.drop(["chr", "start", "end", "ensembl"], axis=1)
        
        total_df = pd.merge(total_df, heterochromatin_merge_df,  how='outer', left_index=True, right_index=True)    
        total_df.to_csv("/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/total_readCount.txt" , sep="\t", header=True,  index=True)
        
        
        
        
        
        
        
        
        
        
        
        #print(total_df[0:5])
        #print(len(total_df))
        
        
        
        
        

if __name__ == "__main__":
    gene_cal = cal_gene_marks()
    #gene_cal.parse_ID()
    #gene_cal.get_gene_body()
    #gene_cal.run_multiBamSummary(markers=["H3K27me3"], bs="H3K27me3.bed", bed_file="/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_body_2k.bed")
    """
    qsub -q batch -V -l nodes=1:ppn=26 ./run_multiBamSummary_H3K27me3.bed.pbs
    """
    #gene_cal.get_TSS(distence=2000)
    #gene_cal.run_multiBamSummary(markers=["H3K4me3"], bs="H3K4me3.bed", bed_file="/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_tss_2000.bed")
    """
    qsub -q batch -V -l nodes=1:ppn=26 ./run_multiBamSummary_H3K4me3.bed.pbs
    """
    #gene_cal.get_TSS(distence=10000)
    #gene_cal.run_multiBamSummary(markers=["H3K9me2", "H3K9me3"], bs="K92_K93.bed", bed_file="/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_tss_10000.bed")
    """
    qsub -q batch -V -l nodes=1:ppn=26 ./run_multiBamSummary_K92_K93.bed.pbs
    """
    #gene_cal.get_TSS(distence=100000)
    #gene_cal.run_multiBamSummary(markers=["H3K27ac", "H3K4me1"], bs="enhancer.bed", bed_file="/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_tss_100000.bed")
    """
    qsub -q batch -V -l nodes=1:ppn=26 ./run_multiBamSummary_enhancer.bed.pbs
    """
    gene_cal.merge_read_count()    
    
    
    