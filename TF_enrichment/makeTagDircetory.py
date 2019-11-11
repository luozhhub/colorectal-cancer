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
    weeks = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
    markers = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
    output_dir = "/home/zhluo/zhao/nfr"
    
    for we in weeks:
        for mark in markers:
            tagDir = os.path.join(tagFiles_dir, "%s_%s" % (we, mark))
            inputDir = os.path.join(tagFiles_dir, "%s_%s" % (we, "Input"))
            output_path = os.path.join(tagDir, "peaks.1kb.txt")
            cmd = "findPeaks %s -i %s -o %s -nfr -size 1000" % (tagDir, inputDir, output_path)
            
            handle = open("/home/zhluo/zhao/pbs/findpeaks/findpeaks_%s_%s.pbs" % (we, mark), "w")
            handle.write(cmd)
            handle.close()


if __name__=="__main__":
    #makeTagDirectory()
    """
    for i in `ls ./*.pbs` ;do qsub -q batch $i ; done
    """
    findPeaks()
    
                
