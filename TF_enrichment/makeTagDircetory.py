"""

"""
import pandas as pd
import os,sys

def makeTagDirectory():
    bamfile_dir = "/home/zhluo/zhao/bamFiles"
    markers = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3", "Input"]
    weeks = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
    replicate = ["1", "2" ,"3"]
    suffix = "_combined_1_paired.merged.nodup.bam"
    output_dir = "/home/zhluo/zhao/create_tags"
    for we in weeks:
        for mark in markers:
            samples = []
            for rep in replicate:
                sample = "%s/%s-%s-%s%s" %(bamfile_dir, we, rep, mark, suffix)
                samples.append(sample)
            output_path = os.path.join(output_dir, "%s_%s" % (we, mark))
            bamfiles = " ".join(samples)
            cmd = "makeTagDirectory %s %s" % (output_path, bamfiles)
            handle = open("/home/zhluo/zhao/pbs/%s_%s.pbs" % (we, mark), "w")
            handle.write(cmd)
            handle.close()           
    
def findPeaks():
    tagFiles_dir = "/home/zhluo/zhao/create_tags"
    weeks = ["2weeks", "4weeks", "7weeks", "10weeks"]
    markers = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
    output_dir = "/home/zhluo/zhao/nfr"
    
    for we in weeks:
        for mark in markers:
            tagDir = os.path.join(tagFiles_dir, "%s_%s" % (we, mark))
            inputDir = os.path.join(tagFiles_dir, "%s_%s" % ("ctrl", mark))
            output_path = os.path.join(tagDir, "peaks.200.txt")
            cmd = "findPeaks %s -i %s -o %s -nfr -size 200" % (tagDir, inputDir, output_path)
            
            handle = open("/home/zhluo/zhao/pbs/findpeaks/findpeaks_%s_%s.pbs" % (we, mark), "w")
            handle.write(cmd)
            handle.close()
            
def add_100bp():
    nfr_dir = "/home/zhluo/zhao/create_tags"
    weeks = ["2weeks", "4weeks", "7weeks", "10weeks"]
    markers = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
    
    for we in weeks:
        for mark in markers:
            peak_path = os.path.join(nfr_dir, "%s_%s" % (we, mark), "peaks.200.txt")
            df = pd.read_csv(peak_path, sep="\t", index_col=False, names=["PeakID", "chr", "start", "end", "strand", "Normalized Tag Count", "Not used",\
                    "findPeaks Score Total Tags", "Control Tags (normalized to IP Experiment)", "Fold Change vs Control", "p-value vs Control" , "Fold Change vs Local", "p-value vs Local" , "Clonal Fold Change"], comment='#')
            df["start"] = [item - 100 for item in df["start"]]
            df["end"] = [item + 100 for item in df["end"]]
            df.drop(["PeakID"], inpulace=True)
            df.to_csv(os.path.join(nfr_dir, "%s_%s_peaks.400.txt" % (we, mark)), sep="\t", header=True,  index=False)
            
            
def transfer_file():
    nfr_dir = "/home/zhluo/zhao/create_tags"
    weeks = ["2weeks", "4weeks", "7weeks", "10weeks"]
    markers = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
    
    for we in weeks:
        for mark in markers:
            mv_files = os.path.join(nfr_dir, "%s_%s_peaks.400.txt" % (we, mark))
            cmd = "scp zhluo@211.69.141.147:%s /home/zhihl/Project/CRC/Chip_analysis/TF_fold_change/nfr_400_peaks/" % (mv_files)
            os.system(cmd)
    


if __name__=="__main__":
    #step1
    #makeTagDirectory()
    """
    for i in `ls ./*.pbs` ;do qsub -q batch $i ; done
    """
    #step2
    #findPeaks()
    #step3
    add_100bp()
    #step4:
    transfer_file()
    
    
                
