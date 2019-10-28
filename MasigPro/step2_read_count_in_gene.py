import sys, os, re
import gzip
import pandas as pd

def add_peakID():
    peak_dir = "/home/zhluo/Project/CRC/data_nazhang/step37_tagAlign/every_sample/"
    files = os.listdir(peak_dir)
    output_dir = "/home/zhluo/Project/CRC/data_nazhang/step37_tagAlign/every_sample_add_peakID/"
    if not os.path.exists(output_dir): 
        os.makedirs(output_dir)
    for one_file in files:
        print(one_file)
        sample_name = one_file.split("_")[0]
        handle = gzip.open(os.path.join(peak_dir, one_file), "rb")
        output = open(os.path.join(output_dir, sample_name + ".bed"), "w")
        n = 0
        for line in handle:
            n = n + 1
            peakID = "peak" + str(n)
            line = line.decode("utf-8").replace("N", peakID)
            output.write(line)
        output.close()
    
def run_annotatePeaks():
    file_dir = "/home/zhluo/Project/CRC/data_nazhang/step37_tagAlign/every_sample_add_peakID"
    output_dir = "/home/zhluo/Project/CRC/data_nazhang/step37_tagAlign/annotatePeaks_output"
    for one_file in os.listdir(file_dir):
        sample_name = one_file.split(".")[0]
        pbs_handle = open(sample_name + ".pbs", "w")
        input_file_path = os.path.join(file_dir, one_file)
        outputfile = os.path.join(output_dir, sample_name + ".anno.bed")
        pbs_handle.write("annotatePeaks.pl %s mm10 > %s" % (input_file_path, outputfile))    
        
        
def run_multiBamSummary(bs = 10000):
    #calculate read count
    #bam files in /home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles
    bamDir = "/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles"
    output_dir = "/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/multiBamSummary_readcount"
    weeks = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
    replicate = ["1", "2" ,"3"]
    markers = ["H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9me2", "H3K9me3"]
    suffix = "_combined_1_paired.merged.nodup.bam"
    sample_list = []    
    lable_list = []
    for we in weeks:
        for mark in markers:
            for rep in replicate:
                sample = "%s-%s-%s%s" %(we, rep, mark, suffix)
                sample_list.append(os.path.join(bamDir, sample))
                lable_list.append("%s_%s_%s" % (we, rep, mark)) 
    bamfiles = " ".join(sample_list)
    lables = " ".join(lable_list)
    cmd = "multiBamSummary bins --binSize %s --numberOfProcessors 25 --bamfiles %s  --minMappingQuality 30 --labels %s -out %s --outRawCounts %s" %(bs, bamfiles, \
    lables, os.path.join(output_dir, "readCounts_%s.npz" % bs), os.path.join(output_dir, "readCounts_%s.tab" % bs))
    pbs_handle = open("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/run_multiBamSummary_%s.pbs" % bs, "w")
    pbs_handle.write(cmd)
    pbs_handle.close()
    
def run_bamCompare():
    bamDir = "/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles"
    weeks = ["2weeks", "4weeks", "7weeks", "10weeks"]
    replicate = ["1", "2" ,"3"]
    markers = ["H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9me2", "H3K9me3"]
    suffix = "_combined_1_paired.merged.nodup.bam"
    sample_list = []    
    lable_list = []
    for we in weeks:
        for mark in markers:
            for rep in replicate:
                sample = "%s-%s-%s%s" %(we, rep, mark, suffix)
                sample_list.append(os.path.join(bamDir, sample))
                
                ctrl_sample =  "%s-%s-%s%s" %(we, rep, mark, suffix)
    
def modify_read_count_file():
    #remove the "#" first
    df = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/multiBamSummary_readcount/readCounts_10000bp.tab", header=0, sep="\t", quotechar="'")
    
    df_sub = df.drop(df.columns[[0,1, 2]], axis=1)
    #print(df_sub[0:5])
    peaks = ["peak%s" % index for index, row in df_sub.iterrows()]
    #print(peaks)
    df_sub.index = peaks
    #df_sub.insert(0, 'peak', peaks)
    print(df_sub[0:5])
    
    df_sub.to_csv("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/multiBamSummary_readcount/expression_read_count_10000.txt", header=True, index=True)

def anno_expression():
    peak_file = "/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/multiBamSummary_readcount/peak_file_10000.txt"
    df = pd.read_csv(peak_file, sep="\t", header=0)
    peaks = ["peak%s" % index for index, row in df.iterrows()]
    df["peak_id"] = peaks
    df["score"] = [1000 for i in range(0, len(df))]
    df["strand"] = ["+" for i in range(0, len(df))]
    print(df[0:5])
    df.to_csv("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/multiBamSummary_readcount/expression_anno_input_10000.txt", sep="\t", header=False, index=False)
    
    pbs_handle = open("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/annotate_output/expression.pbs", "w")
    input_file_path = "/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/multiBamSummary_readcount/expression_anno_input_10000.txt"
    outputfile = os.path.join("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/annotate_output", "expression_10000.anno.bed")
    pbs_handle.write("annotatePeaks.pl %s mm10 > %s" % (input_file_path, outputfile)) 
    pbs_handle.close()

def parse_gene_read_count():
    annofile = "/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/annotate_output/expression_10000.anno.bed"
    df = pd.read_csv(annofile, sep= "\t", header=0)
    peakID = df.iloc[: , 0]
    peak_gene = df["Nearest Ensembl"]
    print(peakID[0:10])
    print(peak_gene[0:10])
    peak_dict = {}
    for i in range(0, len(peakID)):
        if peak_gene[i] is "NaN":
            continue
        peak_dict[peakID[i]] = peak_gene[i]    
    df_readCount = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/multiBamSummary_readcount/expression_read_count_10000.txt", header=0, index_col=0, sep=",")
    gene_name_list  = []
    for index, row in df_readCount.iterrows():
        if index in peak_dict:
            gene_name = peak_dict[index]
            gene_name_list.append(gene_name)
        else: 
            gene_name = "NaN"    
            gene_name_list.append(gene_name)
    df_readCount.index = gene_name_list
    print(df_readCount[0:15])
    df_readCount.to_csv("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/annotate_output/tmp_ensembl.txt", header=True, index=True)
    df_readCount = df_readCount.groupby(df_readCount.index).sum()
    print(df_readCount[0:15])
    
    df_readCount.to_csv("/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/annotate_output/expression_10000_ensembl.txt", header=True, index=True)

if __name__ == "__main__":
    #add_peakID()
    """
    cd /home/zhluo/Project/CRC/data_nazhang/step37_tagAlign/every_sample_add_peakID
    for i in `ls ./` ;do annotatePeaks.pl $i mm10 > "../annotatePeaks_output/"${i%.*}".txt" ;done
    """
    #run_annotatePeaks() 
    """
    cd /home/zhluo/Project/CRC/data_nazhang/step37_tagAlign
    for i in `ls ./*.pbs` ;do qsub $i -l mem=30G;done
    """  
    #run_multiBamSummary(bs=1000)
    """
    qsub -q batch -V -l nodes=1:ppn=26 ./run_multiBamSummary_1000.pbs
    awk -F"\t" '{print $1"\t"$2"\t"$3}' readCounts_10000bp.tab >peak_file_10000.txt
    """
    #modify_read_count_file()
    #anno_expression()
    """
    qsub -l mem=40G ./expression.pbs
    """
    parse_gene_read_count()    