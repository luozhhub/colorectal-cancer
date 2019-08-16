"""
This script used to select chip-seq data from cromwell directory.

"""
#!/usr/bin/python
import os
import re


#tranfer psudo replicate bed file, this is used to create master file
def transfer_psudo_replicate():
    logDir = "/public/home/zhluo/project/CRC_data/step24_cromwell/last_version"
    dataBaseDir = "/public/home/zhluo/project/CRC_data/step24_cromwell/cromwell-executions/chip/"
    files = os.listdir(logDir)
    for oneF in files:
        #get job ID
        if not re.search(".out", oneF):
            continue
        handle = open(oneF, "r")
        lines = handle.readlines()
        jobID = (lines[195].split("/")[9])
    
        #get ppr1 path
        ppr1_dir = os.path.join(dataBaseDir, jobID, "call-pool_ta_pr1/execution")
        files_p = os.listdir(ppr1_dir)
        for oneP in files_p:
            if not re.search(".gz", oneP):
                continue
            bedF = oneP
        ppr1_bed = os.path.join(ppr1_dir, bedF)
        
        #get ppr2 path
        ppr2_dir = os.path.join(dataBaseDir, jobID, "call-pool_ta_pr2/execution")
        files_p = os.listdir(ppr2_dir)
        for oneP in files_p:
            if not re.search(".gz", oneP):
                continue
            bedF = oneP
        ppr2_bed = os.path.join(ppr2_dir, bedF)
    
        #log
        class_t = bedF.split("-")[0]       
        mark = bedF.split("-")[2].split("_")[0]
        print("%s\t%s\t%s" % (class_t, mark, ppr1_bed))
        print("%s\t%s\t%s" % (class_t, mark, ppr2_bed))
        
        #transfer
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step26_pseudo_bed/bedfiles/" % ppr1_bed)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step26_pseudo_bed/bedfiles/" % ppr2_bed)


#transfer bigwig file in pooled sample result 
def transfer_bigWig():
    logDir = "/public/home/zhluo/project/CRC_data/step24_cromwell/last_version"
    dataBaseDir = "/public/home/zhluo/project/CRC_data/step24_cromwell/cromwell-executions/chip/"
    files = os.listdir(logDir)
    for oneF in files:
        #get job ID
        if not re.search(".out", oneF):
            continue
        handle = open(os.path.join(logDir, oneF), "r")
        lines = handle.readlines()
        jobID = (lines[195].split("/")[9])

        #get bigwig path    
        pooled_dir = os.path.join(dataBaseDir, jobID, "call-macs2_pooled/execution")
        files_p = os.listdir(pooled_dir)
        for oneP in files_p:
            if not re.search(".fc.signal.bigwig", oneP):
                continue
            bigwigF = oneP
        pooled_bigwig = os.path.join(pooled_dir, bigwigF)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step28_checkbigWig/bigwig/" % pooled_bigwig)


#transfer pooled bed files, which contained all enriched peaks
def transfer_peak():
    logDir = "/public/home/zhluo/project/CRC_data/step24_cromwell/last_version"
    dataBaseDir = "/public/home/zhluo/project/CRC_data/step24_cromwell/cromwell-executions/chip/"
    files = os.listdir(logDir)
    for oneF in files:
        #get job ID
        if not re.search(".out", oneF):
            continue
        handle = open(os.path.join(logDir, oneF), "r")
        lines = handle.readlines()
        jobID = (lines[195].split("/")[9])

        #get bigwig path    
        pooled_dir = os.path.join(dataBaseDir, jobID, "call-macs2_pooled/execution")
        files_p = os.listdir(pooled_dir)
        for oneP in files_p:
            if not re.search(".narrowPeak.gz", oneP):
                continue
            peakF = oneP
        pooled_peak = os.path.join(pooled_dir, peakF)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step30_pooled_peak/" % pooled_peak)

if __name__ == "__main__":
    transfer_peak()
