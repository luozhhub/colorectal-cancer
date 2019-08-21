"""
first transfer file
scp zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step32_element/enhancer/enhancer.peak.unique.ID.bed 
then run mapping
"""


#!/usr/bin/python
import pandas as pd
import os,sys
import numpy as np
import copy

def transfer_file():
    for ele in ["enhancer", "promoter", "repressed", "heterochromatin"]:
      cmd = "scp zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step32_element/%s/%s.peak.unique.ID.bed /home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/%s/" %(ele, ele, ele)
      os.system(cmd)
    


def mapping_gene_name_to_peak():
    mm_gff = "/home/zhihl/Project/CRC/Chip_analysis/mm10.Tss.gff"
    peak_file_dir = "/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821"
    gene_anno = pd.read_table(mm_gff, names=["chr", "tss", "strand", "gene_type", "gene_id", "trans_type", "trans_id"])
    element = ["enhancer", "promoter", "repressed", "heterochromatin"]
    for ele in element:
        peak_file = pd.read_csv("%s/%s/%s.peak.unique.ID.bed" % (peak_file_dir, ele, ele), names=["chr", "start", "end", "ID"], sep="\t")
        output_file = "%s/%s/%s.peak.unique.ID.gene_name.bed" % (peak_file_dir, ele, ele)
        peak_file["middle"] = (peak_file["end"] - peak_file["start"])/2 + peak_file["start"]
        new_df = pd.DataFrame(columns=peak_file.columns)
        for index,row in peak_file.iterrows():
            try:
                chrom = row["chr"]
                middle = int(row["middle"])
                gene_chr_df = copy.deepcopy(gene_anno[gene_anno["chr"] == chrom])
                gene_chr_df["middle"] = middle
                exit(1)
                gene_chr_df["devi"] = gene_chr_df["middle"] - gene_chr_df["tss"]
                gene_chr_df["devi"] = np.abs(gene_chr_df["devi"])
                idx = gene_chr_df["devi"].idxmin()
                gene_id = gene_chr.loc[idx,"gene_id"]
                row["gene_id"] = gene_id.split(".")[0]
                new_df = new_df.append(row)
            except:
                continue
        new_df.to_csv(output_file, sep="\t", index=0)
        
if __name__ == "__main__":
    #transfer_file()
    mapping_gene_name_to_peak()
    """   
    dirname = "/home/zhihl/Project/CRC/Chip_analysis"
    filename_1 = "mm10.Tss.gff"
    filename_2 = sys.argv[1]
    gene_anno = pd.read_table(os.path.join(dirname, filename_1), header=None)
    gene_anno.columns = ["chr", "tss", "strand", "gene_type", "gene_id", "trans_type", "trans_id"]
    peak_file = pd.read_table(os.path.join(dirname, filename_2), names=["chr", "start", "end"])
    peak_file["middle"] = (peak_file["end"] - peak_file["start"])/2 + peak_file["start"]
    cols = peak_file.columns
    #print (cols)
    #cols.append("gene_id")
    new_df = pd.DataFrame(columns=cols)
    for index,row in peak_file.iterrows():
        try:
            chrom = row["chr"]
            middle = int(row["middle"])
            gene_chr = gene_anno[gene_anno["chr"] == chrom]
            #print (gene_chr[0:5])
            #print (type(middle))
            gene_chr["middle"] = middle
            gene_chr["devi"] = gene_chr["middle"] - gene_chr["tss"]
            gene_chr["devi"] = np.abs(gene_chr["devi"])
            #print (gene_chr)
            idx = gene_chr["devi"].idxmin()
            gene_id = gene_chr.loc[idx,"gene_id"]
            row["gene_id"] = gene_id.split(".")[0]
            new_df = new_df.append(row)
            #print("%s\t%s\t%s\t%s"%(chrom, row["start"], row["end"], row["gene_id"]))
        except:
            continue
    new_df.to_csv(sys.argv[2], sep="\t", index=0)
    """
    