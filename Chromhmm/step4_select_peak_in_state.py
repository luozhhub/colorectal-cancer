from step3_create_master_list import chromhmm_analysis
import os

class element_peaks(chromhmm_analysis):
    def __init__(self):
      chromhmm_analysis.__init__(self)
        
    def select_element_peaks(self):
        element_dir = "/home/zhluo/Project/CRC/data_nazhang/step32_element/"
        master_list_dir = "/home/zhluo/Project/CRC/data_nazhang/step31_master_list"
        
        #enhancer
        H3K27ac = os.path.join(master_list_dir, "H3K27ac.bed.sort.merged")
        enhancer_region = os.path.join(element_dir, "enhancer/enhancer.state.bed.sort.merged")
        enhancer_peak = os.path.join(element_dir, "enhancer/enhancer.peak.bed")
        enhancer_unique_peak = os.path.join(element_dir, "enhancer/enhancer.peak.unique.bed")
        cmd = "bedtools intersect -a %s -b %s -wa  > %s" \
        %(H3K27ac, enhancer_region, enhancer_peak)
        self.run(cmd)
        cmd = "sort %s |uniq > %s" %(enhancer_peak, enhancer_unique_peak)
        self.run(cmd)
        
        #promoter
        H3K4me3 = os.path.join(master_list_dir, "H3K4me3.bed.sort.merged")
        promoter_region = os.path.join(element_dir, "promoter/promoter.state.bed.sort.merged")
        promoter_peak = os.path.join(element_dir, "promoter/promoter.peak.bed")
        promoter_unique_peak = os.path.join(element_dir, "promoter/promoter.peak.unique.bed")
        cmd = "bedtools intersect -a %s -b %s -wa  > %s" \
        %(H3K4me3, promoter_region, promoter_peak)
        self.run(cmd)
        cmd = "sort %s |uniq > %s" %(promoter_peak, promoter_unique_peak)
        self.run(cmd)
        
        #repressed
        H3K27me3 = os.path.join(master_list_dir, "H3K27me3.bed.sort.merged")
        repressed_region = os.path.join(element_dir, "repressed/repressed.state.bed.sort.merged")
        repressed_peak = os.path.join(element_dir, "repressed/repressed.peak.bed")
        repressed_unique_peak = os.path.join(element_dir, "repressed/repressed.peak.unique.bed")
        cmd = "bedtools intersect -a %s -b %s -wa  > %s" \
        %(H3K27me3, repressed_region, repressed_peak)
        self.run(cmd)
        cmd = "sort %s |uniq > %s" %(repressed_peak, repressed_unique_peak)
        self.run(cmd)
        
        #heterochromatin
        H3K9me3 = os.path.join(master_list_dir, "H3K9me3.bed.sort.merged")
        heterochromatin_region = os.path.join(element_dir, "heterochromatin/heterochromatin.state.bed.sort.merged")
        heterochromatin_peak = os.path.join(element_dir, "heterochromatin/heterochromatin.peak.bed")
        heterochromatin_unique_peak = os.path.join(element_dir, "heterochromatin/heterochromatin.peak.unique.bed")
        cmd = "bedtools intersect -a %s -b %s -wa  > %s" \
        %(H3K27me3, heterochromatin_region, heterochromatin_peak)
        self.run(cmd)
        cmd = "sort %s |uniq > %s" %(heterochromatin_peak, heterochromatin_unique_peak)
        self.run(cmd)
        
        
        
    
if __name__ == "__main__":
    elementPeaks = element_peaks()
    elementPeaks.select_element_peaks()