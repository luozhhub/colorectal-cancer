"""
1. transfer tagAlign file to 147, the files in call-bam2ta
      transfer_tagAlign()  

"""


def transfer_tagAlign():
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

        
        
        
        #get every peak path
        every_dir = os.path.join(dataBaseDir, jobID, "call-bam2ta/shard-0/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("tagAlign.gz", oneP):
                continue
            peakF = oneP
        every_0_peak = os.path.join(every_dir, peakF)
        
        every_dir = os.path.join(dataBaseDir, jobID, "call-bam2ta/shard-1/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("tagAlign.gz", oneP):
                continue
            peakF = oneP
        every_1_peak = os.path.join(every_dir, peakF)
        
        every_dir = os.path.join(dataBaseDir, jobID, "call-bam2ta/shard-2/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("tagAlign.gz", oneP):
                continue
            peakF = oneP
        every_2_peak = os.path.join(every_dir, peakF)

        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step37_tagAlign/every_sample" % every_0_peak)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step37_tagAlign/every_sample" % every_1_peak)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step37_tagAlign/every_sample" % every_2_peak)   
        
def transfer_bamFile():
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

        
        
        
        #get every peak path
        every_dir = os.path.join(dataBaseDir, jobID, "call-filter/shard-0/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("nodup.bam$", oneP):
                continue
            peakF = oneP
        every_0_peak = os.path.join(every_dir, peakF)
        print(every_0_peak)
        
        every_dir = os.path.join(dataBaseDir, jobID, "call-filter/shard-1/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("nodup.bam$", oneP):
                continue
            peakF = oneP
        every_1_peak = os.path.join(every_dir, peakF)
        print(every_1_peak)
        
        every_dir = os.path.join(dataBaseDir, jobID, "call-filter/shard-2/execution/")
        files_p = os.listdir(every_dir)
        for oneP in files_p:
            if not re.search("nodup.bam$", oneP):
                continue
            peakF = oneP
        every_2_peak = os.path.join(every_dir, peakF)
        print(every_2_peak)

        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles" % every_0_peak)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles" % every_1_peak)
        os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step38_every_bam/bamFiles" % every_2_peak) 
        
if __name__ == "__main__":
    #transfer_tagAlign()
    transfer_bamFile()
