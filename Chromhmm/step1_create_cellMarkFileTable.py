"""
this script is used to create cellMarkFileTable
usage:
python ./step1_create_cellMarkFileTable.py
cd /home/zhluo/Project/CRC/data_nazhang/pbs/
qsub -q batch -V -l nodes=1:ppn=1 -l mem=60G /home/zhluo/Project/CRC/data_nazhang/pbs/run_chromhmm.pbs
"""
#!/usr/bin/python

import sys,os

#try1: all files
def create_cellMarkFileTable(skip_ctrl=False, outputFile=None):
    dirPath = "/home/zhluo/Project/CRC/data_nazhang/step11_bed"
    file_list = os.listdir(dirPath)
    output= open(outputFile, "w")
    for file in file_list:
        array = file.split("-")
        cellType = array[0]
        mark = array[2].split(".")[0]
        if array[2].split(".")[-1] == "bai":
            continue
        if skip_ctrl is True:
            if cellType == "ctrl":
                continue
        control = cellType + "-" + array[1] + "-" + "Inpu.rmdup.unique.bed"
        if mark == "Inpu":
            continue
        else:
            #output.write("%s\t%s\t%s\n" % (cellType, mark, os.path.join(dirPath, file)))
            output.write("%s\t%s\t%s\t%s\n" % (cellType, mark, file, control))
    output.close()
    
    
    out = open("/home/zhluo/Project/CRC/data_nazhang/pbs/run_chromhmm.pbs", "w")
    binary_cmd = "java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /home/nazhang/luozhihui/software/ChromHMM/CHROMSIZES/mm10.txt /home/zhluo/Project/CRC/data_nazhang/step11_bed /home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/Chromhmm/cellMarkFile_all_sample.txt /home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/binary_files"
    learnModel = " && java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar LearnModel -p 25 /home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/binary_files  /home/zhluo/Project/CRC/data_nazhang/step29_chromhmm/state13 13 mm10"
    out.write(binary_cmd)
    out.write(learnModel)
    out.close()


###############################################################################################
#try every  week data
#that will be no meaning
###############################################################################################
def calculate_each_week():
    #try2weeks
    dirPath = "/home/zhluo/Project/CRC/data_nazhang/step11_bed"
    file_list = os.listdir(dirPath)
    output= open("/home/zhluo/Project/CRC/data_nazhang/cellmarkfiletable_2weeks.txt", "w")
    for file in file_list:
        array = file.split("-")
        cellType = array[0]
        mark = array[2].split(".")[0]
        if array[2].split(".")[-1] == "bai":
            continue
        if cellType != "2weeks":
            continue
        if array[1] == "3" and mark == "H3K9me3":
            continue
        control = cellType + "-" + array[1] + "-" + "Inpu.rmdup.unique.bed"
        if mark == "Inpu":
            continue
        else:
            #output.write("%s\t%s\t%s\n" % (cellType, mark, os.path.join(dirPath, file)))
            output.write("%s\t%s\t%s\t%s\n" % (cellType, mark, file, control))
    output.close()
    binary_cmd = "java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /home/nazhang/luozhihui/software/ChromHMM/CHROMSIZES/mm10.txt /home/zhluo/Project/CRC/data_nazhang/step11_bed /home/zhluo/Project/CRC/data_nazhang/cellmarkfiletable_2weeks.txt /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/2weeks"
    learnModel = " && java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar LearnModel -p 25 /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/2weeks  /home/zhluo/Project/CRC/data_nazhang/step13_chromhmmOutput/2weeks 15 mm10"
    out = open("/home/zhluo/Project/CRC/data_nazhang/pbs/binarizebed_2weeks.pbs", "w")
    out.write(binary_cmd)
    out.write(learnModel)
    out.close()
    
    
    
    #try4weeks
    dirPath = "/home/zhluo/Project/CRC/data_nazhang/step11_bed"
    file_list = os.listdir(dirPath)
    output= open("cellmarkfiletable_4weeks.txt", "w")
    for file in file_list:
        array = file.split("-")
        cellType = array[0]
        mark = array[2].split(".")[0]
        if array[2].split(".")[-1] == "bai":
            continue
        if cellType != "4weeks":
            continue
        control = cellType + "-" + array[1] + "-" + "Inpu.rmdup.unique.bed"
        if mark == "Inpu":
            continue
        else:
            #output.write("%s\t%s\t%s\n" % (cellType, mark, os.path.join(dirPath, file)))
            output.write("%s\t%s\t%s\t%s\n" % (cellType, mark, file, control))
    output.close()
    binary_cmd = "java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /home/nazhang/luozhihui/software/ChromHMM/CHROMSIZES/mm10.txt /home/zhluo/Project/CRC/data_nazhang/step11_bed /home/zhluo/Project/CRC/data_nazhang/cellmarkfiletable_4weeks.txt /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/4weeks"
    learnModel = " && java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar LearnModel -p 25 /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/4weeks  /home/zhluo/Project/CRC/data_nazhang/step13_chromhmmOutput/4weeks 15 mm10"
    out = open("/home/zhluo/Project/CRC/data_nazhang/pbs/binarizebed_4weeks.pbs", "w")
    out.write(binary_cmd)
    out.write(learnModel)
    out.close()
    
    
    #try7weeks
    dirPath = "/home/zhluo/Project/CRC/data_nazhang/step11_bed"
    file_list = os.listdir(dirPath)
    output= open("cellmarkfiletable_7weeks.txt", "w")
    for file in file_list:
        array = file.split("-")
        cellType = array[0]
        mark = array[2].split(".")[0]
        if array[2].split(".")[-1] == "bai":
            continue
        if cellType != "7weeks":
            continue
        control = cellType + "-" + array[1] + "-" + "Inpu.rmdup.unique.bed"
        if mark == "Inpu":
            continue
        else:
            #output.write("%s\t%s\t%s\n" % (cellType, mark, os.path.join(dirPath, file)))
            output.write("%s\t%s\t%s\t%s\n" % (cellType, mark, file, control))
    output.close()
    binary_cmd = "java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /home/nazhang/luozhihui/software/ChromHMM/CHROMSIZES/mm10.txt /home/zhluo/Project/CRC/data_nazhang/step11_bed /home/zhluo/Project/CRC/data_nazhang/cellmarkfiletable_7weeks.txt /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/7weeks"
    learnModel = " && java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar LearnModel -p 25 /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/7weeks  /home/zhluo/Project/CRC/data_nazhang/step13_chromhmmOutput/7weeks 15 mm10"
    out = open("/home/zhluo/Project/CRC/data_nazhang/pbs/binarizebed_7weeks.pbs", "w")
    out.write(binary_cmd)
    out.write(learnModel)
    out.close()
    
    #try10weeks
    dirPath = "/home/zhluo/Project/CRC/data_nazhang/step11_bed"
    file_list = os.listdir(dirPath)
    output= open("cellmarkfiletable_10weeks.txt", "w")
    for file in file_list:
        array = file.split("-")
        cellType = array[0]
        mark = array[2].split(".")[0]
        if array[2].split(".")[-1] == "bai":
            continue
        if cellType != "10weeks":
            continue
        
        if cellType == "10weeks" and array[1] == "2" and mark == "H3K27me3":
            continue
        if cellType == "10weeks" and array[1] == "2" and mark == "H3K4me3":
            continue
        if cellType == "10weeks" and array[1] == "1" and mark == "H3K9me2":
            continue
        control = cellType + "-" + array[1] + "-" + "Inpu.rmdup.unique.bed"
        if mark == "Inpu":
            continue
        else:
            #output.write("%s\t%s\t%s\n" % (cellType, mark, os.path.join(dirPath, file)))
            output.write("%s\t%s\t%s\t%s\n" % (cellType, mark, file, control))
    output.close()
    binary_cmd = "java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /home/nazhang/luozhihui/software/ChromHMM/CHROMSIZES/mm10.txt /home/zhluo/Project/CRC/data_nazhang/step11_bed /home/zhluo/Project/CRC/data_nazhang/cellmarkfiletable_10weeks.txt /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/10weeks"
    learnModel = " && java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar LearnModel -p 25 /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/10weeks  /home/zhluo/Project/CRC/data_nazhang/step13_chromhmmOutput/10weeks 15 mm10"
    out = open("/home/zhluo/Project/CRC/data_nazhang/pbs/binarizebed_10weeks.pbs", "w")
    out.write(binary_cmd)
    out.write(learnModel)
    out.close()
    
    #ctrl
    dirPath = "/home/zhluo/Project/CRC/data_nazhang/step11_bed"
    file_list = os.listdir(dirPath)
    output= open("cellmarkfiletable_ctrl.txt", "w")
    for file in file_list:
        array = file.split("-")
        cellType = array[0]
        mark = array[2].split(".")[0]
        if array[2].split(".")[-1] == "bai":
            continue
        if cellType != "ctrl":
            continue
        
        if cellType == "ctrl" and array[1] == "1" and mark == "H3K9me3":
            continue
        control = cellType + "-" + array[1] + "-" + "Inpu.rmdup.unique.bed"
        if mark == "Inpu":
            continue
        else:
            #output.write("%s\t%s\t%s\n" % (cellType, mark, os.path.join(dirPath, file)))
            output.write("%s\t%s\t%s\t%s\n" % (cellType, mark, file, control))
    output.close()
    binary_cmd = "java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar BinarizeBed -b 200 /home/nazhang/luozhihui/software/ChromHMM/CHROMSIZES/mm10.txt /home/zhluo/Project/CRC/data_nazhang/step11_bed /home/zhluo/Project/CRC/data_nazhang/cellmarkfiletable_ctrl.txt /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/ctrl"
    learnModel = " && java -Xms60000M -jar /home/nazhang/luozhihui/software/ChromHMM/ChromHMM.jar LearnModel -p 25 /home/zhluo/Project/CRC/data_nazhang/step12_chromhmm/ctrl  /home/zhluo/Project/CRC/data_nazhang/step13_chromhmmOutput/ctrl 15 mm10"
    out = open("/home/zhluo/Project/CRC/data_nazhang/pbs/binarizebed_ctrl.pbs", "w")
    out.write(binary_cmd)
    out.write(learnModel)
    out.close()



if __name__ == "__main__":
    create_cellMarkFileTable(outputFile="cellMarkFile_all_sample.txt")