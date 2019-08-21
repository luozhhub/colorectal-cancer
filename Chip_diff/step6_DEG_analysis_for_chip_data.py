"""
use the script "merge_DEG_table.R"
firstly, we should change the count table format, because the dataframe of python and R are not compatibility
just delete the first column in first line is OK

usage:
python step6_DEG_analysis_for_chip_data.py
"""
import os,sys
from subprocess import *

class DEG_chip():
    def __init__(self):
        self.merge_DEG_table = "/home/zhihl/Project/CRC/rna_seq_git/Chip_diff/merge_DEG_table.R"
    
    
    def run(self, cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p
        
    def run_DEG(self):
        marker = ["H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"]
        count_table_dir = "/home/zhihl/Project/CRC/Chip_analysis/dff/version0821"
        for one_mark in marker:
            count_table = os.path.join(count_table_dir, one_mark + "_count_table.txt")
            out_table_name = "all_diff_data_%s.txt" % one_mark
            cmd  = "Rscript %s %s %s %s" %(self.merge_DEG_table, count_table, count_table_dir, out_table_name)
            self.run(cmd)



if __name__ == "__main__":
    DEGChip = DEG_chip()
    DEGChip.run_DEG()            
            

