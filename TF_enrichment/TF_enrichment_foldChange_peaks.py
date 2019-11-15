import os,sys


def obtain_nfr_region():
    nfr_dir = "/home/zhihl/Project/CRC/Chip_analysis/TF_fold_change/nfr_400_peaks"
    weeks = ["2weeks", "4weeks", "7weeks", "10weeks"]
    markers = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
    element = ["enhancer", "enhancer", "promoter", "repressed", "heterochromatin", "heterochromatin"]
    #markers = ["H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"]
    #element = ["enhancer", "promoter", "repressed",  "heterochromatin"]
    element_dict = dict(zip(markers, element))
    for we in weeks:
        for mark in markers:
            nfr_file = os.path.join(nfr_dir, "%s_%s_peaks.400.txt" % (we, mark))
            fold_change_file = os.path.join("/home/zhihl/Project/CRC/Chip_analysis/TF_fold_change", "%s_%s_%s.bed" % (element_dict[mark], mark, we))
            output_file = "/home/zhihl/Project/CRC/Chip_analysis/TF_fold_change/fold_change_nfr/%s_%s.400.bed" % (we, mark)
            cmd = "bedtools intersect -a %s -b %s -wa > %s" % (nfr_file, fold_change_file, output_file)
            os.system(cmd)
            
            
def transfer():
    weeks = ["2weeks", "4weeks", "7weeks", "10weeks"]
    markers = ["H3K27ac", "H3K4me1", "H3K4me3", "H3K27me3", "H3K9me2", "H3K9me3"]
    element = ["enhancer", "enhancer", "promoter", "repressed", "heterochromatin", "heterochromatin"]
    #markers = ["H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"]
    #element = ["enhancer", "promoter", "repressed",  "heterochromatin"]
    element_dict = dict(zip(markers, element))
    for we in weeks:
        for mark in markers:
            peak_file = "/home/zhihl/Project/CRC/Chip_analysis/TF_fold_change/fold_change_nfr/%s_%s.400.bed" % (we, mark)
            cmd = "scp %s zhluo@211.69.141.147:/home/zhluo/zhao/TF_finding/" %(peak_file)
            os.system(cmd)
    
    
    


if __name__ == "__main__":
    #obtain_nfr_region()
    transfer()