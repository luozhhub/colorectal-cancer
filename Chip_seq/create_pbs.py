#!/usr/bin/python


weeks = ["2weeks", "4weeks", "7weeks", "10weeks", "ctrl"]
histon = ["H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9me2", "H3K9me3"]
control = "Input"

for wek in weeks:
    for his in histon:
        json_file = "%s_%s.json" % (wek, his)
        metadata_file = "metadata_%s_%s.json" % (wek, his)
        cmd = "source activate encode-chip-seq-pipeline;cd /public/home/zhluo/project/CRC_data/step24_cromwell; java -jar -Dconfig.file=/public/home/zhluo/project/CRC/chip-seq-pipeline2/backends/backend.conf -Dbackend.providers.Local.config.concurrent-job-limit=2 /public/home/zhluo/cromwell-34.jar run /public/home/zhluo/project/CRC/chip-seq-pipeline2/chip_bwa.wdl -i /public/home/zhluo/project/CRC_data/step24_cromwell/input_json/%s -m /public/home/zhluo/project/CRC_data/step24_cromwell/metadata/%s\n" % (json_file, metadata_file)
        fileout = open("pbs/%s_%s.lsf" % (wek, his), "w")
        header = """#BSUB -J chipJob
#BSUB -q normal
##BSUB -R "span[ptile=20]"
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -o stdout_%J.out
##BSUB -e stderr_%J.err

"""
        fileout.write(header)
        fileout.write(cmd)
        fileout.close()
