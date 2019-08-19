"""
this script is used to select enhancer and promoter from chromhmm result, and then used to create a master list
useage:
#crate merged state
cd /home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/Chromhmm
python ./step3_create_master_list.py

note: something si wrong with the ppr selected peak file
3 works, this script do:
1. process function, merge state
2. state_length_summary function, summary the region length
3. create_master_list function, create the master list
"""

#!/usr/bin/python
import os
import pandas as pd
import sys
from subprocess import *
import re

class chromhmm_analysis():
    def __init__(self):
        self.state_dir = "/home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/state13"
        self.peakfiles_dir = "/home/zhluo/Project/CRC/data_nazhang/step17_macs2_every"
        self.select_peaks_for_every_sample = "/home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/select_peaks_for_every_sample"
        self.select_peaks_for_psudo_sample = "/home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/select_peaks_for_psudo_sample"
        
        self.heterochromatin_states = ["stateE1.bed.sort.merged", "stateE3.bed.sort.merged"]
        self.repressed_states = ["stateE4.bed.sort.merged", "stateE5.bed.sort.merged"]
        self.enhancer_states = ["stateE7.bed.sort.merged", "stateE8.bed.sort.merged", "stateE9.bed.sort.merged"]
        self.promoter_states = ["stateE6.bed.sort.merged", "stateE10.bed.sort.merged", "stateE11.bed.sort.merged", "stateE12.bed.sort.merged", "stateE13.bed.sort.merged"]
        self.state_list = ["E1", "E2", "E3", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E11", "E12", "E13"]
        
        #preprocess
        self.select_peaks_for_every_sample_list = []
        for onefile in os.listdir("/home/zhluo/Project/CRC/data_nazhang/step30_peak/every_peak"):
            if re.search("narrowPeak.gz", onefile):
                self.select_peaks_for_every_sample_list.append(os.path.join("/home/zhluo/Project/CRC/data_nazhang/step30_peak/every_peak", onefile))
                
        self.select_peaks_for_psudo_sample_list = []
        for onefile in os.listdir("/home/zhluo/Project/CRC/data_nazhang/step30_peak/psudo_peak"):
            if re.search("narrowPeak.gz", onefile):
                self.select_peaks_for_psudo_sample_list.append(os.path.join("/home/zhluo/Project/CRC/data_nazhang/step30_peak/psudo_peak", onefile))           

    def run(self, cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p


    def curate_state(self, state=None, outputFile=None ):
        """
        this function is used to create the total region of each state
        """
        bedFiles = ["10weeks_13_segments.bed", "7weeks_13_segments.bed", "4weeks_13_segments.bed", "2weeks_13_segments.bed"]
        df_total  = pd.DataFrame(columns = ["chr", "start", "end", "state"])
        for oneFile in bedFiles:
           file_path = os.path.join(self.state_dir, oneFile)
           df = pd.read_csv(file_path, sep="\t", header=None)
           df.columns = ["chr", "start", "end", "state"]
           sub_df = df.loc[df["state"] == state, : ]
           df_total = df_total.append(sub_df)
        df_total.drop(['state'], axis=1)
        df_total.to_csv(outputFile, sep="\t", header=False, float_format=None, index=False)


    def sort_bed(self, inputFile=None, outputFile=None):
        cmd = "sort -k1,1 -k2,2n %s > %s" %(inputFile, outputFile)
        self.run(cmd=cmd)


    def merge_bed(self, inputFile=None, outputFile=None):
        cmd = "bedtools merge -i %s >%s" %(inputFile, outputFile)
        self.run(cmd=cmd)


    def select_peaks_in_state(self, peakFile=None, state_bed=None, select_peak_files=None):
        #cmd= "bedtools intersect -a /home/zhluo/Project/CRC/data_nazhang/step17_macs2_every/2weeks-2-H3K27ac_peaks.narrowPeak -b /home/zhluo/Project/CRC/data_nazhang/step13_chromhmmOutput/state8.sort.merge.bed -wa > 2weeks-2-H3K27ac_peaks.narrowPeak.state8Anno.bed"
        cmd= "bedtools intersect -a %s -b %s -wa > %s" % (peakFile, state_bed, select_peak_files)
        self.run(cmd=cmd)

    def merge_all_peaks(self, select_peak_files=None, state_total_file=None, zcat=True):
        files_str = " ".join(select_peak_files)
        if zcat is True:
            cmd = "zcat %s > %s" % (files_str, state_total_file)
        else:
            cmd = "cat %s > %s" % (files_str, state_total_file)
        self.run(cmd)
        
        sort_file = state_total_file + ".sort"
        self.sort_bed(inputFile = state_total_file, outputFile = sort_file)

        merged_file = sort_file + ".merged"
        self.merge_bed(inputFile = sort_file, outputFile = merged_file)

    def process(self, state_outDir=None): 
        """
        this function is used to create state file
        """
        #1. merge state in each state
        
        for one_state in self.state_list:
            #create each state bed file
            state_combined_path = "%s/state%s.bed"% (state_outDir, one_state)
            chromhmm_obj.curate_state(state=one_state, outputFile=state_combined_path)
            #sort
            state_combined_sort_file = state_combined_path + ".sort"
            chromhmm_obj.sort_bed(inputFile=state_combined_path, outputFile = state_combined_sort_file)
            #merge
            state_combined_merge_file = state_combined_sort_file + ".merged"
            chromhmm_obj.merge_bed(inputFile = state_combined_sort_file, outputFile = state_combined_merge_file)
        
            
        #2. select peaks 
        ##in each sample
        stateFiles = os.listdir(state_outDir)
        df_state_peaks = pd.DataFrame()
        for stateF in stateFiles:
            if not re.search(".merged", stateF):
                continue
            for peakF in self.select_peaks_for_every_sample_list:
                
                filename = os.path.basename(peakF)
                selectpeakF = os.path.join(self.select_peaks_for_every_sample,  filename.split("_")[0] + "." + stateF.replace(".sort.merged", ""))
                state_F = os.path.join(state_outDir, stateF)
                chromhmm_obj.select_peaks_in_state(peakFile=peakF, state_bed=state_F, select_peak_files=selectpeakF)
                count = len(open(selectpeakF, 'r').readlines())
                df_state_peaks.loc[filename.split("_")[0] , stateF.replace(".sort.merged", "")] = count
        df_state_peaks.to_csv("/home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/select_peaks_for_every_sample/state_peak.txt", sep="\t")
        
        ## in psudo sample
        df_state_peaks = pd.DataFrame()
        for stateF in stateFiles:
            if not re.search(".merged", stateF):
                continue
            for peakF in self.select_peaks_for_psudo_sample_list:
                filename = os.path.basename(peakF)
                selectpeakF = os.path.join(self.select_peaks_for_psudo_sample,  filename.split("_")[0] + "_" + filename.split("_")[3]  + "." + stateF.replace(".sort.merged", ""))
                state_F = os.path.join(state_outDir, stateF)
                chromhmm_obj.select_peaks_in_state(peakFile=peakF, state_bed=state_F, select_peak_files=selectpeakF)
                count = len(open(selectpeakF, 'r').readlines())
                df_state_peaks.loc[filename.split("_")[0] , stateF.replace(".sort.merged", "")] = count
        df_state_peaks.to_csv("/home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/select_peaks_for_psudo_sample/state_peak.txt", sep="\t")
        
    def state_length_summary(self, state_outDir=None):
        """   
        #this function is used to calculate state length in each state
        """    
        state_files = os.listdir(state_outDir)
        for stateF in state_files:
            if not re.search(".merged", stateF):
                continue
            stateName = stateF.split(".")[0]
            df_state = pd.read_csv(os.path.join(state_outDir, stateF), sep="\t", header=None)
            df_state.columns = ["chr", "start", "end"]
            df_state["length"] = df_state["end"] - df_state["start"]
            total_length = sum(df_state["length"])
            print(stateName, float(total_length)/2700000000)
            
            
    def create_master_list(self):
        """
        this function used to create master list
        """    
        pooled_peak_dir = "/home/zhluo/Project/CRC/data_nazhang/step30_peak/pooled_peak"
        peak_list = []
        for one_f in os.listdir(pooled_peak_dir):
            if not re.search("narrowPeak.gz", one_f):
                continue
            peak_list.append(one_f)
        
        histon = ["H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9me2", "H3K9me3"]
        for his in histon:
            select_peak_file = []
            for peak_f in peak_list:
                if re.search(his, peak_f):
                    select_peak_file.append(os.path.join(pooled_peak_dir, peak_f))
            #print(select_peak_file)
            outFile = "/home/zhluo/Project/CRC/data_nazhang/step31_master_list/" + his +".bed"
            chromhmm_obj.merge_all_peaks(select_peak_files = select_peak_file, state_total_file = outFile)
            
        stateDir = "/home/zhluo/Project/CRC/data_nazhang/step31_master_list/"
        state_files = os.listdir(stateDir)
        for stateF in state_files:
            if not re.search(".merged", stateF):
                continue
            stateName = stateF.split(".")[0]
            df_state = pd.read_csv(os.path.join(stateDir, stateF), sep="\t", header=None)
            df_state.columns = ["chr", "start", "end"]
            df_state["length"] = df_state["end"] - df_state["start"]
            total_length = sum(df_state["length"])
            print(stateName, float(total_length)/2700000000)
            
    def element_peak(self, state_outDir=None, element_dir = None):
        stateDir = state_outDir
        #for enhancer
        enhancer_files = [os.path.join(stateDir, mark_s) for mark_s in self.enhancer_states]
        out_enhancer = os.path.join(element_dir, "enhancer/enhancer.state.bed")
        if not os.path.exists(os.path.dirname(out_enhancer)):
            os.makedirs(os.path.dirname(out_enhancer))
        chromhmm_obj.merge_all_peaks(select_peak_files = enhancer_files, state_total_file = out_enhancer, zcat=False)
        
        #for promoter
        promoter_files = [os.path.join(stateDir, mark_s) for mark_s in self.promoter_states]
        out_promoter = os.path.join(element_dir, "promoter/promoter.state.bed")
        if not os.path.exists(os.path.dirname(out_promoter)):
            os.makedirs(os.path.dirname(out_promoter))
        chromhmm_obj.merge_all_peaks(select_peak_files = promoter_files, state_total_file = out_promoter, zcat=False)
        
        #for heterochromatin
        heterochromatin_files = [os.path.join(stateDir, mark_s) for mark_s in self.heterochromatin_states]
        out_heterochromatin = os.path.join(element_dir, "heterochromatin/heterochromatin.state.bed")
        if not os.path.exists(os.path.dirname(out_heterochromatin)):
            os.makedirs(os.path.dirname(out_heterochromatin))
        chromhmm_obj.merge_all_peaks(select_peak_files = heterochromatin_files, state_total_file = out_heterochromatin, zcat=False)
        
        #for repressed
        repressed_files = [os.path.join(stateDir, mark_s) for mark_s in self.repressed_states]
        out_repressed = os.path.join(element_dir, "repressed/repressed.state.bed")
        if not os.path.exists(os.path.dirname(out_repressed)):
            os.makedirs(os.path.dirname(out_repressed))
        chromhmm_obj.merge_all_peaks(select_peak_files = repressed_files, state_total_file = out_repressed, zcat=False)

class primary_statistic():

    def __init__(self, merged_file):
        self.merged_file = merged_file

    def simply_result(self, state_merge_file):
        #df = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/step13_chromhmmOutput/state10.sort.merge.bed", sep="\t", header=None)
        df = pd.read_csv(state_merge_file, sep="\t", header=None)
        df.columns = ["chr", "start", "end"]
        df["length"] = df["end"] - df["start"]
        total_length = sum(df["length"])
        print(total_length, float(total_length)/2700000000)



if __name__ == "__main__":
    #step 1: create whole file
    step = 3
    chromhmm_obj = chromhmm_analysis()    
    #chromhmm_obj.process(state_outDir="/home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/state_merge_region")
    #chromhmm_obj.state_length_summary(state_outDir="/home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/state_merge_region")
    #chromhmm_obj.create_master_list()
    chromhmm_obj.element_peak(state_outDir="/home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/state_merge_region", element_dir="/home/zhluo/Project/CRC/data_nazhang/step32_element")
    """

    
    #create enhancer and promoter state
    if step < 4:
        stateDir = "/home/zhluo/Project/CRC/data_nazhang/step22_state_peaks/state"

        #enhancer_states = ["stateE7.bed.sort.merged", "stateE8.bed.sort.merged", "stateE10.bed.sort.merged"]
        #promoter_states = ["stateE11.bed.sort.merged", "stateE12.bed.sort.merged", "stateE13.bed.sort.merged"]
        enhancer_files = [os.path.join(stateDir, mark_s) for mark_s in self.enhancer_states]
        promoter_files = [os.path.join(stateDir, mark_s) for mark_s in self.promoter_states]
        
        out_enhancer = "/home/zhluo/Project/CRC/data_nazhang/step22_state_peaks/enhancer/enhancer.state.bed"
        chromhmm_obj.merge_all_peaks(select_peak_files = enhancer_files, state_total_file = out_enhancer)
        out_promoter = "/home/zhluo/Project/CRC/data_nazhang/step22_state_peaks/promoter/promoter.state.bed" 
        chromhmm_obj.merge_all_peaks(select_peak_files = promoter_files, state_total_file = out_promoter)
        
        
    if step < 5 :
    """
