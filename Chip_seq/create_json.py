#!/usr/bin/python


weeks = ["2weeks", "4weeks", "7weeks", "10weeks", "ctrl"]
histon = ["H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9me2", "H3K9me3"]
control = "Input"

for wek in weeks:
    for his in histon:
        file_content="""
        {
        "use_bwa_mem_for_pe" : true,
        "chip.pipeline_type" : "histone",
        "chip.genome_tsv" : "/public/home/zhluo/project/CRC_data/step24_cromwell/genome_database/mm10.tsv",
        "chip.fastqs" : [
            [["/public/home/zhluo/project/CRC_data/chip_gz/%s-1-%s_combined_1_paired.fq.gz",
              "/public/home/zhluo/project/CRC_data/chip_gz/%s-1-%s_combined_2_paired.fq.gz"]],
            [["/public/home/zhluo/project/CRC_data/chip_gz/%s-2-%s_combined_1_paired.fq.gz",
              "/public/home/zhluo/project/CRC_data/chip_gz/%s-2-%s_combined_2_paired.fq.gz"]],
            [["/public/home/zhluo/project/CRC_data/chip_gz/%s-3-%s_combined_1_paired.fq.gz",
              "/public/home/zhluo/project/CRC_data/chip_gz/%s-3-%s_combined_2_paired.fq.gz"]]
         ],
        "chip.ctl_fastqs" : [
            [["/public/home/zhluo/project/CRC_data/chip_gz/%s-1-Input_combined_1_paired.fq.gz",
              "/public/home/zhluo/project/CRC_data/chip_gz/%s-1-Input_combined_2_paired.fq.gz"]],
            [["/public/home/zhluo/project/CRC_data/chip_gz/%s-2-Input_combined_1_paired.fq.gz",
              "/public/home/zhluo/project/CRC_data/chip_gz/%s-2-Input_combined_2_paired.fq.gz"]],
            [["/public/home/zhluo/project/CRC_data/chip_gz/%s-3-Input_combined_1_paired.fq.gz",
              "/public/home/zhluo/project/CRC_data/chip_gz/%s-3-Input_combined_2_paired.fq.gz"]]
        ],

        "chip.paired_end" : true,

        "chip.always_use_pooled_ctl" : true,
        "chip.title" : "histon modification",
        "chip.description" : "histon ChIP-seq on CRC",

        "chip.bwa_cpu" : 10,
        "chip.bwa_mem_mb" : 20000,
        "chip.bwa_time_hr" : 24,

        "chip.filter_cpu" : 10,
        "chip.filter_mem_mb" : 16000,
        "chip.filter_time_hr" : 24,

        "chip.bam2ta_cpu" : 10,
        "chip.bam2ta_mem_mb" : 16000,
        "chip.bam2ta_time_hr" : 24,

        "chip.spr_mem_mb" : 16000,

        "chip.fingerprint_cpu" : 10,
        "chip.fingerprint_mem_mb" : 16000,
        "chip.fingerprint_time_hr" : 24,

        "chip.xcor_cpu" : 10,
        "chip.xcor_mem_mb" : 16000,
        "chip.xcor_time_hr" : 24,

        "chip.macs2_mem_mb" : 20000,
        "chip.macs2_time_hr" : 24,

        "chip.spp_cpu" : 10,
        "chip.spp_mem_mb" : 20000,
        "chip.spp_time_hr" : 72
        }
        """ % (wek, his, wek, his, wek, his, wek, his, wek, his, wek, his, wek, wek, wek, wek, wek, wek)
        outfile = open("input_json/%s_%s.json" % (wek, his), "w")
        outfile.write(file_content)
        outfile.close()




