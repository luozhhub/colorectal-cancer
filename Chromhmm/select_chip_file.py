"""
this script is used to transfer cromwell out file to 147 cluster.
"""

#!/usr/bin/python
import os
import re


#transfer psudoreplicate file
logDir = "/public/home/zhluo/project/CRC_data/step24_cromwell"
dataBaseDir = "/public/home/zhluo/project/CRC_data/step24_cromwell/cromwell-executions/chip/"
files = os.listdir(logDir)
for oneF in files:
    if not re.search(".out", oneF):
        continue
    #print(oneF)
    handle = open(oneF, "r")
    lines = handle.readlines()
    jobID = (lines[195].split("/")[9])

    ppr1_dir = os.path.join(dataBaseDir, jobID, "call-pool_ta_pr1/execution")
    files_p = os.listdir(ppr1_dir)
    for oneP in files_p:
        if not re.search(".gz", oneP):
            continue
        bedF = oneP
    ppr1_bed = os.path.join(ppr1_dir, bedF)
    
    ppr2_dir = os.path.join(dataBaseDir, jobID, "call-pool_ta_pr2/execution")
    files_p = os.listdir(ppr2_dir)
    for oneP in files_p:
        if not re.search(".gz", oneP):
            continue
        bedF = oneP
    ppr2_bed = os.path.join(ppr2_dir, bedF)

    #print(ppr1_bed)
    class_t = bedF.split("-")[0]
    
    mark = bedF.split("-")[2].split("_")[0]
    print("%s\t%s\t%s" % (class_t, mark, ppr1_bed))
    print("%s\t%s\t%s" % (class_t, mark, ppr2_bed))
    #os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step26_pseudo_bed/bedfiles/" % ppr1_bed)
    os.system("scp %s zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step26_pseudo_bed/bedfiles/" % ppr2_bed)

    
