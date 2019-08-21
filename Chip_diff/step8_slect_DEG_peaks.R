#--------------install package-----------------------------------
install_biomaRt = function (){
  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("biomaRt", version = "3.8")
}

#--------------analysis-------------------------------------------
library("biomaRt")
library("org.Hs.eg.db")

#chip_data
###this script used to output DEG peaks, so we need to manually change df_all and setwd for 4 element dir
df_all = read.delim("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/all_diff_data_H3K9me3.txt",sep = "\t", header=T)
setwd("/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/heterochromatin/")
#old chip_daa
##df_all = read.delim("/home/zhihl/Project/CRC/Chip_analysis/dff/all_diff_data_promoter.txt",sep = "\t", header=T)
##setwd("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/")
#for RNA
##df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis_0418_non_unique_bam/all_diff_data.txt",sep = "\t", header=T)
mouse_gene_transcript = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_mouse_gene_transcript.txt", sep="\t", header=T)
orthologs_table = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_humanGene_mouseGene.txt", sep="\t", header=T)


#3. convert essemble to entrez 
essemble_to_entrez = function (deg_list){
  entrez_map = as.data.frame(unlist(as.list(org.Mm.egENSEMBL2EG)))
  ind<-match(deg_list, rownames(entrez_map))
  deg_eg = entrez_map[ind,]
  return(deg_eg)
}

#4. convert essemble transcript to gene
trans2gene = function(deg_list){
  mouse_table = mouse_gene_transcript
  rownames(mouse_table) = mouse_table$Transcript.stable.ID
  select_table = mouse_table[deg_list,]
  deg_gene = select_table[,1]
  return(deg_gene)
}


#5. orthologs of mouse and human
mouse2human = function(deg_list){
  orthologs = orthologs_table
  ind = match(deg_list, orthologs$Mouse.gene.stable.ID)
  human_genes = orthologs[ind,]
  gene_list = unique(na.omit(human_genes[,1]))
  return(gene_list)
}

#6. diff orthologs
diff_gene = function(col_1, col_2, fold_change=1 ,id_type="SYMBOL"){
  #col_1: col name of log2 fold change
  #col_2: col name of padj
  #fold_change: log2 fold change
  #id_type: "SYMBOL", "ENTRIZEID" ....
  #return list of up and down genes
  week_df = df_all[,c(col_1, col_2)]
  diff_week_up = unique(rownames(na.omit(week_df[week_df[,1] > fold_change & week_df[,2] < 0.01, ])))
  diff_week_down = unique(rownames(na.omit(week_df[week_df[,1] < -(fold_change) & week_df[,2] < 0.01, ])))
  return (list(diff_week_up, diff_week_down))
}

diff_genes_10w = diff_gene(col_1 = "log2FoldChange.10wkVSctrl", col_2 = "padj.10wkVSctrl")
diff_genes_7w = diff_gene(col_1 = "log2FoldChange.7wkVSctrl", col_2 = "padj.7wkVSctrl")
diff_genes_4w = diff_gene(col_1 = "log2FoldChange.4wkVSctrl", col_2 = "padj.4wkVSctrl")
diff_genes_2w = diff_gene(col_1 = "log2FoldChange.2wkVSctrl", col_2 = "padj.2wkVSctrl")

# divide groups
commom_cancer_up = intersect(diff_genes_10w[[1]], diff_genes_7w[[1]])
commom_cancer_down = intersect(diff_genes_10w[[2]], diff_genes_7w[[2]])

commom_colits_up = intersect(diff_genes_4w[[1]], diff_genes_2w[[1]])
commom_colits_down = intersect(diff_genes_4w[[2]], diff_genes_2w[[2]])

commom_all_up = intersect(commom_cancer_up, commom_colits_up)
commom_all_down = intersect(commom_cancer_down, commom_colits_down)

commom_can_up = setdiff(commom_cancer_up,commom_all_up)
commom_can_down = setdiff(commom_cancer_down, commom_all_down)

commom_col_up = setdiff(commom_colits_up,commom_all_up)
commom_col_down = setdiff(commom_colits_down, commom_all_down)

commom_can_list = list(commom_cancer_up, commom_cancer_down)
commom_col_list = list(commom_colits_up, commom_colits_down)
commom_all_list = list(commom_all_up, commom_all_down)

#write table
write.table(c(commom_can_list[[1]], commom_can_list[[2]]) , "commom_can_list_up.txt", sep="\n", quote=F, row.names=F, col.names = F, na="NA", eol="\n")
write.table(c(commom_col_list[[1]], commom_col_list[[2]]), "commom_col_list_up.txt", sep="\n", quote=F, row.names=F, col.names = F, na="NA", eol="\n")


print(c(length(diff_genes_10w[[1]]), length(diff_genes_10w[[2]])))
print(c(length(diff_genes_7w[[1]]), length(diff_genes_7w[[2]])))
print(c(length(diff_genes_4w[[1]]), length(diff_genes_4w[[2]])))
print(c(length(diff_genes_2w[[1]]), length(diff_genes_2w[[2]])))
print(c(length(commom_cancer_up), length(commom_cancer_down)))
print(c(length(commom_colits_up), length(commom_colits_down)))
print(c(length(commom_all_up), length(commom_all_down)))


