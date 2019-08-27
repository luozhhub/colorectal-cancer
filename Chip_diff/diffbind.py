#!/usr/bin/python
import os
import pandas as pd
import re


def everyfile():
    #pd.Dataframe(columns=["SampleID","Tissue","Factor","Condition","Treatment","Replicate","bamReads","ControlID","bamControl","Peaks","PeakCaller"])
    #df = pd.DataFrame(columns=["SampleID","Tissue","Factor","Condition","Treatment","Replicate","bamReads","bamControl","Peaks","PeakCaller"])
    df = pd.DataFrame(columns=["SampleID","bamReads","bamControl","Peaks"])
    cellTable = "cellmarkfiletable_all.txt"
    df_cell = pd.read_table("/home/zhihl/Project/CRC/cellmarkfiletable_all.txt", header=None)
    df_cell.columns= columns=["cellType", "marker", "inbed", "controlbed"]
    bamdir = "/home/zhihl/Project/CRC/step9_unique"
    peakdir = "/home/zhihl/Project/CRC/step17_macs2_every"
    for index, row in df_cell.iterrows():
        cellType = row["cellType"]
        if cellType == "7weeks" or cellType == "10weeks":
            tissue = "tumor"
        elif cellType == "2weeks" or cellType == "4weeks":
            tissue = "inflam"
        else:
            tissue = "normal"
        factor = "ER"
        condition = cellType
        treatment = row["marker"]
        replicate = row["inbed"].split(".")[0].split("-")[1]
        bamreads = os.path.join(bamdir, row["inbed"].replace(".bed", ".sort.bam"))
        bamcontrol = os.path.join(bamdir, row["controlbed"].replace(".bed", ".sort.bam"))
        Peaks = os.path.join(peakdir, row["inbed"].split(".")[0] + "_peaks.narrowPeak")
        caller = "macs2"
        sers = pd.Series([row["inbed"].split(".")[0],  bamreads, bamcontrol, Peaks],index=["SampleID","bamReads","bamControl","Peaks"])
        df = df.append(sers, ignore_index=True)
    df.to_csv('diffbind_all.txt',index=False,sep=',')



def mergefile():
    files = os.listdir("/home/zhluo/Project/CRC/data_nazhang/step16_macs2")
    bamdir = "/home/zhluo/Project/CRC/data_nazhang/step15_merge"
    bamlist = os.listdir(bamdir)

    for file in bamlist:
        if not re.search("sort.bam", file):
            continue
        if re.research("Inpu.sort.bam", file):
            continue
        sample = file.replace(".sort.bam", "")
        sampleid = sample
        if sample.split("_")[0] == "7weeks" or sample.split("_")[0] == "10weeks":
            tissue = "tumor"
        else:
            tissue = "normal"
        factor = sample.split("_")[1]
        bamreads = file
        bamcontrol = sample.split("_")[0] + "_Inpu.sort.bam"
        
           
    

if __name__ == "__main__":
    everyfile()
