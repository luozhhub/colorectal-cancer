#!/usr/bin/python3

import os
import sys
from subprocess import *

class quality_control():

    def __init__(self):
        self.fastq_dir = "/home/nazhang/luozhihui/project/CRC/Chip/ChIP-seq/"
        self.Trimmomatic = "/home/nazhang/luozhihui/software/Trimmomatic-0.38/trimmomatic-0.38.jar"
        self.outputDir = "/home/nazhang/luozhihui/project/CRC/trimmoResult"
        self.adaptor = "/home/nazhang/luozhihui/project/CRC/TruSeq2-PE.fa"
        self.ref = "/home/zhluo/Project/CRC/data_nazhang/refernece_genome/genecode/GRCm38.primary_assembly.genome.fa"

    def run(self, cmd=None, wkdir=None):
        sys.stderr.write("Running %s ...\n" % cmd)
        p = Popen(cmd, shell=True, cwd=wkdir)
        p.wait()
        return p

    def fastp(self, fastq1=None, fastq2=None, outputDir=None):
        sample = os.path.basename(fastq1).split(".")[0]
        output1 = sample + ".fastp.fq.gz"
        output1 = os.path.join(outputDir, output1)
        sample = os.path.basename(fastq2).split(".")[0]
        output2 = sample + ".fastp.fq.gz"
        output2 = os.path.join(outputDir, output2)
        sample = sample.split("_")[0]
        htmlPath = os.path.join(outputDir, sample + ".html")
        jsonPath = os.path.join(outputDir, sample + ".json")
        cmd = "fastp -z 4 -i %s -I %s -o %s -O %s -h %s -j %s" % (fastq1, fastq2, output1, output2, htmlPath, jsonPath)
        return cmd

    def extract(self, zipFile=None):
        """
        unzip the fastq file

        """

    def trimAdapterByTrimmomatic(self, fastq1=None, fastq2=None, outputDir=None):
        """
        fastq eg. "823-RNA_L8_1.fq.gz"
        ILLUMINACLIP:%s:2:30:10:
        %s is the adaptor file
        2 is bigest mismatch
        30 is palindrome method threshold value
        10表示simple方法的匹配阈值
        LEADING:3 表示切掉reads 5’端(the start of read)质量低于3的碱基或N
        TRAILING:3 表示切掉reads 3’端(the end of read)质量低于3的碱基或N
        SLIDINGWINDOW:4:15 表示以4个碱基作为窗口，窗口一个碱基一个碱基往后移，如果窗口内碱基的平均质量小于15,后面的都切掉
        MINLEN:36 #以上步骤处理后，如果reads的长度小于36，这条reads也会被排除

        :param fastq1: string
        :param fastq2: string
        :return:
        """
        fastq1_name = os.path.basename(fastq1)
        fastq2_name = os.path.basename(fastq2)
        forward_paired = os.path.join(outputDir, fastq1_name.replace("_R1.fastp.fq.gz", "_1_paired.fq.gz"))
        forward_unpaired = os.path.join(outputDir, fastq1_name.replace("_R1.fastp.fq.gz", "_1_unpaired.fq.gz"))
        reverse_paired = os.path.join(outputDir, fastq2_name.replace("_R2.fastp.fq.gz", "_2_paired.fq.gz"))
        reverse_unpaired = os.path.join(outputDir, fastq2_name.replace("_R2.fastp.fq.gz", "_2_unpaired.fq.gz"))
        logFile = os.path.join(outputDir, fastq1_name.replace("_R1.fastp.fq.gz", ".log"))
        cmd = "java -jar %s PE -threads 4 -phred33 -trimlog %s\
         %s  %s %s %s %s %s \
         ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50" % \
              (self.Trimmomatic, logFile, fastq1, fastq2, forward_paired, forward_unpaired, reverse_paired, reverse_unpaired, self.adaptor)
        return cmd

    def fastqc(self, fastq=None, outputDir=None):
        #dirName = os.path.dirname(fastq)
        cmd = "fastqc -t 6 -o %s %s" % (outputDir, fastq )
        return cmd

    def star(self, fastq1=None, fastq2=None, outputDir=None):
        fastq1_name = os.path.basename(fastq1)
        prefix = os.path.join(outputDir, fastq1_name.replace("_1_paired.fq.gz", ""))
        cmd = "STAR --genomeDir star_genome --readFilesIn %s %s \
        --readFilesCommand zcat --outSAMstrandField intronMotif --runThreadN 8 --outFileNamePrefix %s" % \
              (fastq1, fastq2, prefix)
        return cmd

    def run_bwa(self, core=6, ref=None, fastq_1=None, fastq_2=None, outsam=None, rg=None):
        cmd = "%s mem -R '@RG\\tID:1\\tPL:ILLUMINA\\tSM:%s' -t %s %s %s %s > %s"%("bwa", rg, core, ref, fastq_1, fastq_2, outsam)
        return cmd

class count_read():

    def __init__(self):
        self.samtools = "/home/zxchen/anaconda3/bin/samtools"
        self.picard = "/home/zhluo/Software/picard.jar"


    def sam_to_bam(self, samFile=None, outputDir=None):
        fileName = os.path.basename(samFile)
        bamFile = os.path.join(outputDir, fileName.replace(".sam", ".bam"))
        cmd = "samtools view -@ 4 -F 4 -bS %s -o %s" % (samFile, bamFile)
        return cmd

    def bamsort(self, core=4, mem="3G", bamfile=None, tmpdir=None ,sort_by_name=False, outputDir=None):
        fileName = os.path.basename(bamfile)
        sortbamfile = os.path.join(outputDir, fileName.replace(".bam", ".sort.bam"))
        if sort_by_name is True:
            cmd="%s sort -n -@ %s -m %s -T %s -o %s %s"%("samtools", core, mem, tmpdir, sortbamfile, bamfile)
        else:
            cmd="%s sort -@ %s -m %s -T %s -o %s %s"%("samtools", core, mem, tmpdir, sortbamfile, bamfile)
        return cmd

    def Markduplicate(self, tmpdir=None, input_bam=None, marked_duplicated_bam=None, marked_dup_metrics_txt=None):
        cmd="java -Djava.io.tmpdir=%s -jar -Xms1024m -Xmx10240m %s MarkDuplicates I=%s O=%s M=%s VALIDATION_STRINGENCY=LENIENT ASSUME_SORT_ORDER=queryname REMOVE_DUPLICATES=true"\
            %(tmpdir, self.picard, input_bam, marked_duplicated_bam, marked_dup_metrics_txt)
        return cmd

    def unique_mapping(self, input_bam=None, output_dir=None):
        fileName = os.path.basename(input_bam)
        smple = fileName.split(".")[0] 
        headerFile = os.path.join(output_dir, sample + ".header.txt")
        unique_sam_t = os.path.join(output_dir, sample + ".unique.tmp.sam")
        unique_sam =  os.path.join(output_dir, sample + ".unique.sam")
        sort_bam =  os.path.join(output_dir, sample + ".rmdup.unique.sort.bam")
        cmd = "samtools view -h -q 1 -F 4 -F 256 %s |grep -v XA:Z | grep -v SA:Z | samtools view -@ 10 -Sb - | samtools sort -@ 10 >%s" %(input_bam, sort_bam)
        #cmd="samtools view -H %s > %s && samtools view %s |grep NH:i:1 >%s && cat %s %s > %s && samtools view -@ 10 -Sb %s|samtools sort -@ 10 >%s" %(input_bam, headerFile, input_bam, unique_sam_t, headerFile, unique_sam_t, unique_sam, unique_sam, sort_bam)
        return cmd


if __name__ == "__main__":
    import re
    qc = quality_control()
    step = 16
    #step 1, generate fastqc file
    if step < 1:
        fastqFiles = os.listdir(qc.fastq_dir)
        for fi in fastqFiles:
            fastq = os.path.join(qc.fastq_dir, fi)
            cmd = qc.fastqc(fastq=fastq, outputDir="/home/nazhang/luozhihui/project/CRC/Chip/original_fastqc")
            qc.run(cmd=cmd)

    #step 2, run fastp
    if step < 2:
        fastqFiles = os.listdir(qc.fastq_dir)
        sample_num = []
        for file in fastqFiles:
            regx = re.compile("(.*)_combined_(.*)\.fastq")
            result = regx.search(file)
            number = result.groups()[0]
            if number not in sample_num:
                sample_num.append(number)
        for sample in sample_num:
            fastq1 = "%s_combined_R1.fastq" % (sample)
            fastq1 = os.path.join(qc.fastq_dir, fastq1)
            fastq2 = "%s_combined_R2.fastq" % (sample)
            fastq2 = os.path.join(qc.fastq_dir, fastq2)
            cmd = qc.fastp(fastq1=fastq1, fastq2=fastq2, outputDir="/home/nazhang/luozhihui/project/CRC/Chip/fastp_step2")
            #print("cd /home/nazhang/luozhihui/project/CRC/Chip;" + cmd)
            qc.run(cmd=cmd)


    #step 3, run trimomatic
    if step <3:
        fastq_dir = "/home/nazhang/luozhihui/project/CRC/Chip/fastp_step2"
        fastqFiles = os.listdir(fastq_dir)
        if len(fastqFiles) == 0:
            exit(1)
        sample_num = []
        for file in fastqFiles:
            if re.search(".json", file):
                continue
            if re.search(".html", file):
                continue
            regx = re.compile("(.*)_combined_(.*)\.fastp.fq.gz")
            result = regx.search(file)
            number = result.groups()[0]
            if number not in sample_num:
                sample_num.append(number)
        #print(sample_num)
        #exit(1)
        for sample in sample_num:
            fastq1 = "%s_combined_R1.fastp.fq.gz" % (sample)
            fastq1 = os.path.join(fastq_dir, fastq1)
            fastq2 = "%s_combined_R2.fastp.fq.gz" % (sample)
            fastq2 = os.path.join(fastq_dir, fastq2)
            cmd = qc.trimAdapterByTrimmomatic(fastq1=fastq1, fastq2=fastq2, outputDir="/home/nazhang/luozhihui/project/CRC/Chip/trimomatic_step3")
            #qc.run(cmd)
            OP = open("/home/nazhang/luozhihui/project/CRC/Chip/pbs/trimomatic_" + sample + ".pbs", "w")
            OP.write("cd /home/nazhang/luozhihui/project/CRC/Chip;" + cmd + "\n")
            OP.close()

    
    #step 4, run fastqc
    if step < 4:
        fastq_dir = "/home/nazhang/luozhihui/project/CRC/Chip/trimomatic_step3"
        fastqFiles = os.listdir(fastq_dir)
        if len(fastqFiles) == 0:
            exit(1)
        for fi in fastqFiles:
            if re.search("fq.gz", fi):
                fastq = os.path.join(fastq_dir, fi)
                cmd = qc.fastqc(fastq=fastq, outputDir="/home/nazhang/luozhihui/project/CRC/Chip/fastqc_step4")
                #qc.run(cmd=cmd)
                OP = open("/home/nazhang/luozhihui/project/CRC/Chip/pbs/fastqc_" + fi + ".pbs", "w")
                OP.write("cd /home/nazhang/luozhihui/project/CRC/Chip;" + cmd + "\n")
                OP.close()

   
    #step 5, run bwa
    if step < 5:
        fastq_dir = "/home/zhluo/Project/CRC/data_nazhang/trimomatic_step3"
        fastqFiles = os.listdir(fastq_dir)
        if len(fastqFiles) == 0:
            exit(1)

        sample_num = []
        for file in fastqFiles:
            if re.search(".log", file):
                continue
            regx = re.compile("(.*)_combined_(.*)_paired.fq.gz")
            result = regx.search(file)
            if not result:
                continue
            number = result.groups()[0]
            if number not in sample_num:
                sample_num.append(number)

        for sample in sample_num:
            fastq1 = "%s_combined_1_paired.fq.gz" % (sample)
            fastq1 = os.path.join(fastq_dir, fastq1)
            fastq2 = "%s_combined_2_paired.fq.gz" % (sample)
            fastq2 = os.path.join(fastq_dir, fastq2)
            sam = os.path.join("/home/zhluo/Project/CRC/data_nazhang/step5_bwa_test", sample + ".sam")
            print (sample)
            cmd = qc.run_bwa(core=6, ref=qc.ref, fastq_1=fastq1, fastq_2=fastq2, outsam=sam, rg=sample)
            #cmd = qc.star(fastq1=fastq1, fastq2=fastq2, outputDir="/home/nazhang/luozhihui/project/CRC/star_step5")
            #qc.run(cmd)
            OP = open("/home/zhluo/Project/CRC/data_nazhang/pbs/bwa_" + sample + ".pbs", "w")
            OP.write("cd /home/zhluo/Project/CRC/data_nazhang;" + cmd + "\n")
            OP.close()


    bamPro = count_read()

    #step 6, convert sam to bam
    if step < 6:
        sam_dir = "/home/zhluo/Project/CRC/data_nazhang/step5_bwa"
        samFiles = os.listdir(sam_dir)
        if len(samFiles) == 0:
            exit(1)

        for file in samFiles:
            if not re.search(".sam", file):
                continue
            samfile = os.path.join(sam_dir, file)
            cmd = bamPro.sam_to_bam(samFile=samfile, outputDir="/home/zhluo/Project/CRC/data_nazhang/step6_bam")
            sample = file.strip("\n").strip(".sam")
            OP = open("/home/zhluo/Project/CRC/data_nazhang/pbs/samtools_" + sample + ".pbs", "w")
            OP.write("cd /home/zhluo/Project/CRC/data_nazhang;" + cmd +"\n")
            OP.close()
            #qc.run(cmd)

  
    #step 7, sort bam
    if step < 7:
        bam_dir = "/home/zhluo/Project/CRC/data_nazhang/step6_bam"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            if not re.search(".bam", file):
                continue
            bamfile = os.path.join(bam_dir, file)
            cmd = bamPro.bamsort(bamfile=bamfile, sort_by_name=True, outputDir="/home/zhluo/Project/CRC/data_nazhang/step7_sort", tmpdir="/home/zhluo/Project/CRC/data_nazhang/tmp")
            sample = file.strip("\n").strip(".bam")
            OP = open("/home/zhluo/Project/CRC/data_nazhang/pbs/sort_" + sample + ".pbs", "w")
            OP.write("cd /home/zhluo/Project/CRC/data_nazhang;" + cmd +"\n")
            OP.close()
            

    #step 8, Markduplicate
    if step < 8:
        bam_dir = "/home/zhluo/Project/CRC/data_nazhang/step7_sort"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            if not re.search(".bam", file):
                continue
            bamfile = os.path.join(bam_dir, file)
            sample = file.strip("\n").strip(".sort.bam")
            output_dir = "/home/zhluo/Project/CRC/data_nazhang/step8_mkdup"
            mkdup_file = os.path.join(output_dir, sample + ".mkdup.bam")
            mkmatrix = os.path.join(output_dir, sample + ".matrix.txt")
            cmd = bamPro.Markduplicate(tmpdir="/home/zhluo/Project/CRC/data_nazhang/tmp", input_bam=bamfile, marked_duplicated_bam=mkdup_file, marked_dup_metrics_txt=mkmatrix)
            OP = open("/home/zhluo/Project/CRC/data_nazhang/pbs/mkdup_" + sample + ".pbs", "w")
            OP.write("cd /home/zhluo/Project/CRC/data_nazhang;" + cmd +"\n")
            OP.close()
    

    #step 9, unique mapping reads
    if step < 9:
        bam_dir = "/home/zhluo/Project/CRC/data_nazhang/step8_mkdup"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            if not re.search(".bam", file):
                continue
            bamfile = os.path.join(bam_dir, file)
            #print(file)
            sample = file.strip("\n").replace(".mkdup.bam", "")
            #print(sample)
            output_dir = "/home/zhluo/Project/CRC/data_nazhang/step9_unique"
            #unique_file = os.path.join(output_dir, sample + ".mkdup.unique.bam")
            cmd = bamPro.unique_mapping(input_bam=bamfile, output_dir=output_dir)
            OP = open("/home/zhluo/Project/CRC/data_nazhang/pbs/unique_" + sample + ".pbs", "w")
            OP.write("cd /home/zhluo/Project/CRC/data_nazhang;" + cmd +"\n")
            OP.close()
  

    #step 15, merge bam
    if step < 15:
        bam_dir = "/home/zhluo/Project/CRC/data_nazhang/step9_unique"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)
        bam_dict = {}
        for file in bamFiles:
            #print (file)
            if re.search(".bai", file):
                continue
            sample = file.split(".")[0]
            if sample == "2weeks-3-H3K9me3" or sample == "10weeks-2-H3K27me3" or sample == "10weeks-2-H3K4me3" or sample == "10weeks-1-H3K9me2" or sample == "ctrl-1-H3K9me3":
                continue
            type_feature = sample.split("-")[0] + "_" + sample.split("-")[2]
            bamfile = os.path.join(bam_dir, file)
            if type_feature in bam_dict:
                bam_dict[type_feature].append(bamfile)
            else:
                bam_dict[type_feature] = [bamfile]
                
        for key in bam_dict.keys():
            merge_bam = os.path.join("/home/zhluo/Project/CRC/data_nazhang/step15_merge" , key + ".bam")
            input_bams = " ".join(bam_dict[key])
            #cmd = "samtools merge -@ 10 %s %s" % (merge_bam, input_bams)
            sortBam = merge_bam.replace(".bam", ".sort.bam")
            cmd = "samtools sort -@ 10 -o %s %s" %(sortBam, merge_bam)
            OP = open("/home/zhluo/Project/CRC/data_nazhang/pbs/merge_" + key + ".pbs", "w")
            OP.write("cd /home/zhluo/Project/CRC/data_nazhang;" + cmd +"\n")
            OP.close()

            
    


    #step 16, macs2
    if step < 16:
        bam_dir = "/home/zhluo/Project/CRC/data_nazhang/step15_merge"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            if not re.search(".sort.bam", file):
                continue
            if re.search("Inpu.sort.bam", file):
                continue
            #bamfile = os.path.join(bam_dir, file)
            #cmd = bamPro.bamsort(bamfile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/Chip/sort_step7", tmpdir="/home/nazhang/luozhihui/project/CRC/Chip/tmp")
            sample = file.replace(".sort.bam", "")
            array = sample.split("_")
            array[-1] = "Inpu.sort.bam"
            ctrol = os.path.join(bam_dir , "_".join(array))
            file1 = os.path.join(bam_dir, file)
            OP = open("/home/zhluo/Project/CRC/data_nazhang/pbs/macs2_" + file + ".pbs", "w")
            cmd = "macs2 callpeak -t %s -c %s -f BAMPE -B --nomodel -g mm -q 0.05 -n %s --outdir %s" % (file1, ctrol, sample, "/home/zhluo/Project/CRC/data_nazhang/step16_macs2")
            OP.write("cd /home/nazhang/luozhihui/project/CRC/Chip;" + cmd +"\n")
            OP.close()
            

   

    #step 17, macs2
    if step < 17:
        bam_dir = "/home/zhluo/Project/CRC/data_nazhang/step9_unique"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            if  re.search("rmdup.unique.sort.bam.bai", file):
                continue
            if re.search("Inpu.rmdup.unique.sort.bam", file):
                continue
            #bamfile = os.path.join(bam_dir, file)
            #cmd = bamPro.bamsort(bamfile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/Chip/sort_step7", tmpdir="/home/nazhang/luozhihui/project/CRC/Chip/tmp")
            sample = file.replace(".rmdup.unique.sort.bam", "")
            array = sample.split("-")
            array[-1] = "Inpu.rmdup.unique.sort.bam"
            ctrol = os.path.join(bam_dir , "-".join(array))
            file1 = os.path.join(bam_dir, file)
            OP = open("/home/zhluo/Project/CRC/data_nazhang/pbs/macs2_" + file + ".pbs", "w")
            cmd = "macs2 callpeak -t %s -c %s -f BAMPE -B --nomodel -g mm -q 0.05 -n %s --outdir %s" % (file1, ctrol, sample, "/home/zhluo/Project/CRC/data_nazhang/step17_macs2_every")
            OP.write("cd /home/nazhang/luozhihui/project/CRC/Chip;" + cmd +"\n")
            OP.close()
            
    exit(1)

    #step 8, htseq
    if step < 8:
        bam_dir = "/home/nazhang/luozhihui/project/CRC/sort_step7"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            bamfile = os.path.join(bam_dir, file)
            cmd = bamPro.htseq_count(sortBamFile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/htseq_step8" )
            qc.run(cmd)

    if step < 9:
        bam_dir = "/home/nazhang/luozhihui/project/CRC/sort_step7"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            bamfile = os.path.join(bam_dir, file)
            cmd = bamPro.htseq_count_1(sortBamFile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/htseq_step9")
            qc.run(cmd)

    if step < 10:
        bam_dir = "/home/nazhang/luozhihui/project/CRC/sort_step7"
        bamFiles = os.listdir(bam_dir)
        if len(bamFiles) == 0:
            exit(1)

        for file in bamFiles:
            bamfile = os.path.join(bam_dir, file)
            cmd = bamPro.htseq_count_2(sortBamFile=bamfile, outputDir="/home/nazhang/luozhihui/project/CRC/htseq_step10")
            qc.run(cmd)
