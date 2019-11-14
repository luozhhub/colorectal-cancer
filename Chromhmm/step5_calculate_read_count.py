"""
usage:
1. use transfer_psudo_replicate function in step2, get bed files
2.
"""

#!/usr/bin/python
import os,sys
import re
from subprocess import *
from step4_select_peak_in_state import element_peaks
import pandas as pd
#import pandas as pd

class get_read_count(element_peaks):
    def __init__(self):
        element_peaks.__init__(self)

    def add_peak_id(self):
        """
        this function use to process the master list file, and then add peak ID fro each peak
        the master list can be :/home/zhluo/Project/CRC/data_nazhang/step25_master_list/enhancer/enhancer.peak.bed
        or  /home/zhluo/Project/CRC/data_nazhang/step32_element/enhancer/enhancer.peak.unique.bed
        """
        element = ["enhancer", "promoter", "repressed", "heterochromatin"]
        element_dir = "/home/zhluo/Project/CRC/data_nazhang/step32_element/"
        for ele in element:
            df = pd.read_csv(os.path.join(element_dir, ele + "/" + ele + ".peak.unique.bed"), names=["chr", "start", "end"], sep="\t")
            ser1 = pd.Series(["peak" + str(i) for i in range(len(df))])
            df["peakID"] = ser1
            #print(df[0:5])
            df.to_csv(os.path.join(element_dir, ele + "/" + ele + ".peak.unique.ID.bed"), sep="\t", header=None, index=None)
        

    def read_count(self, element=None, marker=None):
        """
        this function is used to calculate the read count
        the add peak id file can be :/home/zhluo/Project/CRC/data_nazhang/step27_pseudo_diff/peak_name_list/enhancer.peak.id.bed
        or /home/zhluo/Project/CRC/data_nazhang/step32_element/enhancer/enhancer.peak.unique.ID.bed
        """
        outputDir = "/home/zhluo/Project/CRC/data_nazhang/step33_pseudo_diff/read_count_bed"
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        pseudo_dir = "/home/zhluo/Project/CRC/data_nazhang/step26_pseudo_bed/bedfiles/"

        #element = ["enhancer", "promoter", "repressed", "heterochromatin"]
        #marker = ["H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"]
        
        element_dict = dict(zip(element, marker))
        element_dir = "/home/zhluo/Project/CRC/data_nazhang/step32_element/"       
        for ele in element:
            peak_file = os.path.join(element_dir, ele + "/" + ele + ".peak.unique.ID.bed")
            bed_files = os.listdir(pseudo_dir)
            for oneF in bed_files:
                if not re.search(element_dict[ele], oneF):
                    continue
                select_file = os.path.join(pseudo_dir, oneF)            
                if re.search("pr1", oneF):
                    outFName = oneF.split("_")[0] + ".pr1.count.bed"
                else:
                    outFName = oneF.split("_")[0] + ".pr2.count.bed"                    
                outFile = os.path.join(outputDir, outFName)
                cmd = "bedtools intersect -a %s -b %s  -c > %s" % (peak_file, select_file, outFile)
                self.run(cmd)
            
    def merge_samples(self, element=None, marker=None):
        """
        this function is used to merge different sample to one kind of marker table
        
        """
        psudo_bed_dir = "/home/zhluo/Project/CRC/data_nazhang/step33_pseudo_diff/read_count_bed"
        group_dict = {}
        #element = ["enhancer", "promoter", "repressed", "heterochromatin"]
        #marker = ["H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"]
        #element_dict = dict(zip(element, marker))
        
        for ele in marker:
            group_dict[ele] = ["ctrl-1-%s.pr1.count.bed"%ele,  "ctrl-1-%s.pr2.count.bed"%ele,  "2weeks-1-%s.pr1.count.bed"%ele,  "2weeks-1-%s.pr2.count.bed"%ele,  "4weeks-1-%s.pr1.count.bed"%ele,
            "4weeks-1-%s.pr2.count.bed"%ele,  "7weeks-1-%s.pr1.count.bed"%ele,  "7weeks-1-%s.pr2.count.bed"%ele,  "10weeks-1-%s.pr1.count.bed"%ele,  "10weeks-1-%s.pr2.count.bed"%ele]
            df_blank = pd.DataFrame()
            for one_bed in group_dict[ele]:
                print(os.path.join(psudo_bed_dir, one_bed))
                df = pd.read_csv(os.path.join(psudo_bed_dir, one_bed), header=None, index_col=3, sep="\t", names=["chr", "start", "end", "count"])
                df_blank[one_bed] = df["count"]
            df_name = ele + "_count_table.txt"
            df_blank.to_csv(os.path.join("/home/zhluo/Project/CRC/data_nazhang/step33_pseudo_diff/", df_name), sep="\t", header=True, index=True)

if __name__ == "__main__":
    getReadCount = get_read_count()
    #getReadCount.add_peak_id()
    #run in 201908
    """
    getReadCount.read_count(element=["enhancer", "promoter", "repressed", "heterochromatin"], marker=["H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"])
    getReadCount.merge_samples(element=["enhancer", "promoter", "repressed", "heterochromatin"], marker=["H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"])
    """
    getReadCount.read_count(element=["enhancer", "heterochromatin"], marker=["H3K4me1", "H3K9me2"])
    getReadCount.merge_samples(element=["enhancer", "heterochromatin"], marker=["H3K4me1", "H3K9me2"])