"""
use the script "merge_DEG_table.R"
firstly, we should change the count table format, because the dataframe of python and R are not compatibility
just delete the first column in first line is OK

"""

class DEG_chip():
    def __init__(self):
        self.merge_DEG_table = 
        pass
    
    def run(self, cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p
        
    def run_DEG():
        marker = ["H3K27ac", "H3K4me3", "H3K27me3", "H3K9me3"]
        count_table_dir = "/home/zhihl/Project/CRC/Chip_analysis/dff/version0821"
        for one_mark in marker:
            count_table = os.path.join(count_table_dir, one_mark + "_count_table.txt")
            cmd  = "Rscript ./merge_DEG_table.R ./enhancer_count_table.txt ./version0821/ all_diff_data_H3K27ac.txt"
            

