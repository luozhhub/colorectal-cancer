"""
this work is in /home/zhluo/Project/CRC/data_nazhang/step40_masigpro_gene_list/
"""
#!/usr/bin/pyhton
import os,sys
import pandas as pd




def run_computeMatrix():
    bw_dir = "/home/zhluo/Project/CRC/data_nazhang/step28_checkbigWig/bigwig"
    time = ["ctrl", "2weeks", "4weeks", "7weeks", "10weeks"]
    
    #for H3K27ac
    marker = "H3K27ac"
    bed_file = "/home/zhluo/Project/CRC/data_nazhang/colorectal-cancer/MasigPro/gene_tss_1.bed"
    df_enhancer_bed = pd.read_csv(bed_file, sep="\t", header=None, index_col=3, names=["chr", "start", "end"])
    outout_region_dir = "/home/zhluo/Project/CRC/data_nazhang/step40_masigpro_gene_list/computeMatrix_region"
    output_matrix_dir = "/home/zhluo/Project/CRC/data_nazhang/step40_masigpro_gene_list/matrix_output"
    deeptools_chip_indensity_dir = "/home/zhluo/Project/CRC/data_nazhang/step40_masigpro_gene_list/deeptools_chip_intensity"
    output_table_dir = "/home/zhluo/Project/CRC/data_nazhang/step40_masigpro_gene_list/output_table"
    
    
    #for H3K27ac cluster 1
    #zhihl@ubuntu:~/Project/CRC/Chip_analysis/MaSigPro/rna$ scp ./* zhluo@211.69.141.147:/home/zhluo/Project/CRC/data_nazhang/step40_masigpro_gene_list/rna_gene_list/
    for marker in ["H3K27ac", "H3K4me1", "H3K4me3", "H3K9me3", "H3K9me2", "H3K27me3"]:
        for cluster in ["cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9"]:
            gene_cluster_1 = pd.read_csv("/home/zhluo/Project/CRC/data_nazhang/step40_masigpro_gene_list/rna_gene_list/rna_gene_list_%s.txt" % (cluster), sep="\t", index_col=0)
            df_cluster_1 = df_enhancer_bed.loc[gene_cluster_1.index]
            cluster_1_bed_path = os.path.join(outout_region_dir, "RNA_%s_region_%s.bed" % (cluster, marker))
            
            
            #df_cluster_1 = pd.to_numeric(df_cluster_1, downcast='signed')
            df_cluster_1 = df_cluster_1.dropna()
            df_cluster_1 = df_cluster_1.astype({'start': 'int32', 'end': 'int32'})
            
            df_cluster_1.to_csv(cluster_1_bed_path , sep="\t", header=False,  index=False)
            matrix_cluster_1 = os.path.join(output_matrix_dir, "RNA_%s_region_%s.matrix.gz" % (cluster, marker))
            deeptools_chip_intensity_cluster_1 = os.path.join(deeptools_chip_indensity_dir, "RNA_%s_region_%s.pdf" % (cluster, marker))
            output_table_cluster_1  = os.path.join(output_table_dir, "RNA_%s_region_%s.table" % (cluster, marker))
            pbs_handle = open(os.path.join("/home/zhluo/Project/CRC/data_nazhang/step40_masigpro_gene_list/pbs", "RNA_%s_region_%s.pbs" % (cluster, marker)), "w")
            
            
            file_list = []
            samplesLabel = []
            for week in time:
                sample_prefix = "%s-1-%s" % (week, marker)
                for one_file in os.listdir(bw_dir):
                    if sample_prefix in one_file:
                        file_list.append(os.path.join(bw_dir, one_file))
                        samplesLabel.append(sample_prefix)
            
            cmd = "computeMatrix reference-point -S %s -R %s  -p 25 --samplesLabel %s --beforeRegionStartLength 10000  --afterRegionStartLength 10000 --skipZeros -o %s &&  plotProfile -m  %s  --outFileName %s  --outFileNameData  %s --perGroup" % (" ".join(file_list), cluster_1_bed_path, " ".join(samplesLabel), matrix_cluster_1, matrix_cluster_1, deeptools_chip_intensity_cluster_1, output_table_cluster_1)
            pbs_handle.write(cmd)
            pbs_handle.close()
    """
    for i in `ls ./*.pbs` ; do qsub -l nodes=1:ppn=26 $i;done
    """
        
        
    


if __name__ == "__main__":
    run_computeMatrix()