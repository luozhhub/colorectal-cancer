#--------------install package-----------------------------------
install_biomaRt = function (){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt", version = "3.8")
}

#--------------analysis-------------------------------------------
library("biomaRt")
library("org.Mm.eg.db")

df_all = read.delim("/home/zhihl/Project/CRC/Chip_analysis/dff/all_diff_data_enhancer.txt",sep = "\t", header=T)
#df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis_0418_non_unique_bam/all_diff_data.txt", sep = "\t", header=T)
mouse_gene_transcript = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_mouse_gene_transcript.txt", sep="\t", header=T)
orthologs_table = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_humanGene_mouseGene.txt", sep="\t", header=T)
peak_gene = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/enhancer.peak_gene.map.txt", sep="\t", header=T)


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


#------------------------------------------------------------
peak2gene = function(peak_list){
  peak_gene_table = peak_gene
  ind = match(peak_list, peak_gene_table$peakID)
  human_genes = peak_gene_table[ind,]
  gene_list = unique(na.omit(human_genes[,1]))
  return(gene_list)
}



#6. diff orthologs
diff_gene = function(col_1, col_2, fold_change=1 ,id_type="ENSEMBL"){
  #col_1: col name of log2 fold change
  #col_2: col name of padj
  #fold_change: log2 fold change
  #id_type: "SYMBOL", "ENTRIZEID" ....
  #return list of up and down genes
  
  week_df = df_all[,c(col_1, col_2)]
  diff_week_up = rownames(na.omit(week_df[week_df[,1] > fold_change & week_df[,2] < 0.01, ]))
  diff_week_down = rownames(na.omit(week_df[week_df[,1] < -(fold_change) & week_df[,2] < 0.01, ]))
  down_gene = peak2gene(diff_week_down)
  up_gene = peak2gene(diff_week_up)
  
  cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
  up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(up_gene), columns=cols, keytype="ENSEMBL")
  up_list = unique(up_symbol[, id_type])
  down_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(down_gene), columns=cols, keytype="ENSEMBL")
  down_list = unique(down_symbol[, id_type])
  return (list(up_list, down_list))
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

commom_col_list_chip = list(commom_col_up, commom_col_down)
commom_can_list_chip = list(commom_can_up, commom_can_down)
commom_all_list_chip = list(commom_all_up, commom_all_down)


#-----------------------------------------------------------------------------
#RNA
df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis/all_diff_data.txt",sep = "\t", header=T)

diff_gene_in_mouse = function(col_1, col_2, fold_change=1 ,id_type="ENSEMBL"){
  #col_1: col name of log2 fold change
  #col_2: col name of padj
  #fold_change: log2 fold change
  #id_type: "SYMBOL", "ENTRIZEID" ....
  #return list of up and down genes
  
  week_df = df_all[,c(col_1, col_2)]
  diff_week_up = rownames(na.omit(week_df[week_df[,1] > fold_change & week_df[,2] < 0.01, ]))
  diff_week_down = rownames(na.omit(week_df[week_df[,1] < -(fold_change) & week_df[,2] < 0.01, ]))
  
  
  cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
  up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(diff_week_up), columns=cols, keytype="ENSEMBL")
  up_list = unique(up_symbol[, id_type])
  down_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(diff_week_down), columns=cols, keytype="ENSEMBL")
  down_list = unique(down_symbol[, id_type])
  return (list(up_list, down_list))
}

diff_genes_10w = diff_gene_in_mouse("log2FoldChange.10wkVSctrl", "padj.10wkVSctrl")
diff_genes_7w = diff_gene_in_mouse("log2FoldChange.7wkVSctrl", "padj.7wkVSctrl")
diff_genes_4w = diff_gene_in_mouse("log2FoldChange.4wkVSctrl", "padj.4wkVSctrl")
diff_genes_2w = diff_gene_in_mouse("log2FoldChange.2wkVSctrl", "padj.2wkVSctrl")

tummor_vs_control = diff_gene_in_mouse("log2FoldChange.tumorVSctrl", "padj.tumorVSctrl")
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

commom_col_list = list(commom_col_up, commom_col_down)
commom_can_list = list(commom_can_up, commom_can_down)
commom_all_list = list(commom_all_up, commom_all_down)




# only cancer gene
commom_ccan_up = intersect(commom_can_list[[1]], commom_can_list_chip[[1]])
pValue = phyper(length(commom_ccan_up), length(commom_can_list[[1]]), 20000 - length(commom_can_list[[1]]), length(commom_can_list_chip[[1]]), lower.tail = F)
print(paste("up", length(commom_ccan_up), length(commom_can_list[[1]]), 20000 - length(commom_can_list[[1]]), length(commom_can_list_chip[[1]]), pValue, sep=","))

commom_ccan_down = intersect(commom_can_list[[2]], commom_can_list_chip[[2]])
pValue = phyper(length(commom_ccan_down), length(commom_can_list[[2]]), 20000 - length(commom_can_list[[2]]), length(commom_can_list_chip[[2]]), lower.tail = F)
print(paste("down", length(commom_ccan_down), length(commom_can_list[[2]]), 20000 - length(commom_can_list[[2]]), length(commom_can_list_chip[[2]]), pValue, sep=","))

#only colites gene

commom = intersect(commom_col_list[[1]], commom_col_list_chip[[1]])
pValue = phyper(length(commom), length(commom_col_list[[1]]), 20000 - length(commom_col_list[[1]]), length(commom_col_list_chip[[1]]), lower.tail = F)
print(paste("up", length(commom), length(commom_col_list[[1]]), 20000 - length(commom_col_list[[1]]), length(commom_col_list_chip[[1]]), pValue, sep=","))

commom = intersect(commom_col_list[[2]], commom_col_list_chip[[2]])
pValue = phyper(length(commom), length(commom_col_list[[2]]), 20000 - length(commom_col_list[[2]]), length(commom_col_list_chip[[2]]), lower.tail = F)
print(paste("down", length(commom), length(commom_col_list[[2]]), 20000 - length(commom_col_list[[2]]), length(commom_col_list_chip[[2]]), pValue, sep=","))

# all share gene
commom = intersect(commom_all_list[[1]], commom_all_list_chip[[1]])
pValue = phyper(length(commom), length(commom_all_list[[1]]), 20000 - length(commom_all_list[[1]]), length(commom_all_list_chip[[1]]), lower.tail = F)
print(paste("up", length(commom), length(commom_all_list[[1]]), 20000 - length(commom_all_list[[1]]), length(commom_all_list_chip[[1]]), pValue, sep=","))

commom = intersect(commom_all_list[[2]], commom_all_list_chip[[2]])
pValue = phyper(length(commom), length(commom_all_list[[2]]), 20000 - length(commom_all_list[[2]]), length(commom_all_list_chip[[2]]), lower.tail = F)
print(paste("down", length(commom), length(commom_all_list[[2]]), 20000 - length(commom_all_list[[2]]), length(commom_all_list_chip[[2]]), pValue, sep=","))

library("clusterProfiler")
setwd("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/")
ego <- enrichGO(gene          = commom_ccan_up,
                keyType       = "SYMBOL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)
pdf("promoter.caner.up.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

ego <- enrichGO(gene          = commom_ccan_down,
                keyType       = "SYMBOL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)
pdf("promoter.caner.down.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


write.table(commom_can_up, "commom_can_up.txt", sep="\n", quote=F, row.names=T, col.names = T, na="NA", eol="\n")
write.table(commom_ccan_up, "commom_ccan_up.txt", sep="\n", quote=F, row.names=T, col.names = T, na="NA", eol="\n")

human_orthologs = mouse2human(commom_ccan_up)
cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
up_symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(human_orthologs), columns=cols, keytype="ENSEMBL")
up_list = unique(up_symbol[, id_type])

human_orthologs = mouse2human(commom_ccan_down)
cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
down_symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(human_orthologs), columns=cols, keytype="ENSEMBL")
down_list = unique(down_symbol[, id_type]) 

CRCgwasList = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/CRCgwasGene.txt",sep = "\t", header=F)
CRCgwasList = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/CRCmappedGenes.txt",sep = "\t", header=F)

CRCgwasList = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(CRCgwasList$V1), columns=cols, keytype="ENSEMBL")
CRCgwasList = unique(CRCgwasList[, "SYMBOL"])
intersect(CRCgwasList, down_list)