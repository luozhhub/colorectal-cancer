"""
This script used to select chip-seq data from cromwell directory.
useage:
cd /public/home/zhluo/project/CRC_data/step24_cromwell
scp xx@xxxx:/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/Chip_seq/create_path.py ./create_path.py
python ./create_path.py
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

        
        #get pooled peak path    
        pooled_dir = os.path.join(dataBaseDir, jobID, "call-macs2_pooled/execution")
        files_p = os.listdir(pooled_dir)
        for oneP in files_p:
            if not re.search("500K.bfilt.narrowPeak.gz", oneP):
                continue
            peakF = oneP
        pooled_peak = os.path.join(pooled_dir, peakF)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step30_peak/pooled_peak" % pooled_peak)
        
        
        
        #get psudo peak path
        psudo_dir = os.path.join(dataBaseDir, jobID, "call-macs2_ppr1/execution/")
        files_p = os.listdir(psudo_dir)
        for oneP in files_p:
            if not re.search("500K.bfilt.narrowPeak.gz", oneP):
                continue
            peakF = oneP
        psudo_1_peak = os.path.join(psudo_dir, peakF)
        
        psudo_dir = os.path.join(dataBaseDir, jobID, "call-macs2_ppr2/execution/")
        files_p = os.listdir(psudo_dir)
        for oneP in files_p:
            if not re.search("500K.bfilt.narrowPeak.gz", oneP):
                continue
            peakF = oneP
        psudo_2_peak = os.path.join(psudo_dir, peakF)

        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step30_peak/psudo_peak" % psudo_1_peak)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step30_peak/psudo_peak" % psudo_2_peak)  
        
        #get every peak path
        every_dir = os.path.join(dataBaseDir, jobID, "call-macs2/shard-0/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("500K.bfilt.narrowPeak.gz", oneP):
                continue
            peakF = oneP
        every_0_peak = os.path.join(every_dir, peakF)
        
        every_dir = os.path.join(dataBaseDir, jobID, "call-macs2/shard-1/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("500K.bfilt.narrowPeak.gz", oneP):
                continue
            peakF = oneP
        every_1_peak = os.path.join(every_dir, peakF)
        
        every_dir = os.path.join(dataBaseDir, jobID, "call-macs2/shard-2/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("500K.bfilt.narrowPeak.gz", oneP):
                continue
            peakF = oneP
        every_2_peak = os.path.join(every_dir, peakF)

        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step30_peak/every_peak" % every_0_peak)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step30_peak/every_peak" % every_1_peak)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step30_peak/every_peak" % every_2_peak)  
        
        
        
def transfer_pooled_tagAlign():
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
        pooled_dir = os.path.join(dataBaseDir, jobID, "call-pool_ta/execution")
        files_p = os.listdir(pooled_dir)
        for oneP in files_p:
            if not re.search(".tagAlign.gz", oneP):
                continue
            bigwigF = oneP
        pooled_bigwig = os.path.join(pooled_dir, bigwigF)
        os.system("scp zhluo@211.69.141.130:%s /home/zhihl/zhaochen/chip/pooled_tagAlign/" % pooled_bigwig)      
             

if __name__ == "__main__":
    #transfer_peak()
    transfer_pooled_tagAlign()
