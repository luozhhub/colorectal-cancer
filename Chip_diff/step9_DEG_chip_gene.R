#--------------install package-----------------------------------
install_biomaRt = function (){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt", version = "3.8")
}


#--------------analysis-------------------------------------------
library("biomaRt")
library("org.Mm.eg.db")

df_all = read.delim("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/all_diff_data_H3K9me3.txt",sep = "\t", header=T)
#df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis_0418_non_unique_bam/all_diff_data.txt", sep = "\t", header=T)
mouse_gene_transcript = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_mouse_gene_transcript.txt", sep="\t", header=T)
orthologs_table = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_humanGene_mouseGene.txt", sep="\t", header=T)
peak_gene = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/heterochromatin/heterochromatin.peak.unique.ID.gene_name.bed", sep="\t", header=T)


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
peak2gene = function(peak_list, peak_gene_path){
  peak_gene = read.delim(peak_gene_path, sep="\t", header=T)
  peak_gene_table = peak_gene
  ind = match(peak_list, peak_gene_table$ID)
  human_genes = peak_gene_table[ind,]
  gene_list = unique(na.omit(human_genes$gene_id))
  return(gene_list)
}



#6. Chip
diff_gene_for_chip = function(df_all, col_1, col_2, fold_change=1 ,id_type="ENSEMBL", peak_gene){
  #col_1: col name of log2 fold change
  #col_2: col name of padj
  #fold_change: log2 fold change
  #id_type: "SYMBOL", "ENTRIZEID" ....
  #return list of up and down genes
  
  week_df = df_all[,c(col_1, col_2)]
  diff_week_up = rownames(na.omit(week_df[week_df[,1] > fold_change & week_df[,2] < 0.01, ]))
  diff_week_down = rownames(na.omit(week_df[week_df[,1] < -(fold_change) & week_df[,2] < 0.01, ]))
  down_gene = peak2gene(diff_week_down, peak_gene)
  up_gene = peak2gene(diff_week_up, peak_gene)
  
  cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
  up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(up_gene), columns=cols, keytype="ENSEMBL")
  up_list = unique(up_symbol[, id_type])
  down_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(down_gene), columns=cols, keytype="ENSEMBL")
  down_list = unique(down_symbol[, id_type])
  return (list(up_list, down_list))
}




#-----------------------------------------------------------------------------
#for RNA
#df_rna = read.delim("/home/zhihl/Project/CRC/RNA_analysis/all_diff_data.txt",sep = "\t", header=T)

diff_gene_in_mouse = function(df_all, col_1, col_2, fold_change=1 ,id_type="ENSEMBL"){
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


Summary_data = function(mark, state, diff_gene_func, data_type){
  summary = list()
  if (data_type == "RNA")
  {
    df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis/all_diff_data.txt",sep = "\t", header=T)
    diff_genes_10w = diff_gene_func(df_all, "log2FoldChange.10wkVSctrl", "padj.10wkVSctrl", 1, "ENSEMBL")
    diff_genes_7w = diff_gene_func(df_all, "log2FoldChange.7wkVSctrl", "padj.7wkVSctrl", 1, "ENSEMBL")
    diff_genes_4w = diff_gene_func(df_all, "log2FoldChange.4wkVSctrl", "padj.4wkVSctrl", 1, "ENSEMBL")
    diff_genes_2w = diff_gene_func(df_all, "log2FoldChange.2wkVSctrl", "padj.2wkVSctrl", 1, "ENSEMBL")
  }
  if (data_type == "Chip")
  {
    df_all = read.delim(paste("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/all_diff_data_", mark, ".txt", sep=""),sep = "\t", header=T)
    #peak gene mapping
    peak_gene = paste("/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/", state, "/", state, ".peak.unique.ID.gene_name.bed", sep="")
    diff_genes_10w = diff_gene_func(df_all, "log2FoldChange.10wkVSctrl", "padj.10wkVSctrl", 1, "ENSEMBL", peak_gene)
    diff_genes_7w = diff_gene_func(df_all, "log2FoldChange.7wkVSctrl", "padj.7wkVSctrl", 1, "ENSEMBL", peak_gene)
    diff_genes_4w = diff_gene_func(df_all, "log2FoldChange.4wkVSctrl", "padj.4wkVSctrl", 1, "ENSEMBL", peak_gene)
    diff_genes_2w = diff_gene_func(df_all, "log2FoldChange.2wkVSctrl", "padj.2wkVSctrl", 1, "ENSEMBL", peak_gene)
  }
  
  
  tummor_vs_control = diff_gene_in_mouse(df_rna, "log2FoldChange.tumorVSctrl", "padj.tumorVSctrl", 1, "ENSEMBL")
  
  #7week and 10week
  commom_cancer_up = intersect(diff_genes_10w[[1]], diff_genes_7w[[1]])
  commom_cancer_down = intersect(diff_genes_10w[[2]], diff_genes_7w[[2]])
  #2week and 4week
  commom_colits_up = intersect(diff_genes_4w[[1]], diff_genes_2w[[1]])
  commom_colits_down = intersect(diff_genes_4w[[2]], diff_genes_2w[[2]])
  #2,4,7,10 week
  commom_all_up = intersect(commom_cancer_up, commom_colits_up)
  commom_all_down = intersect(commom_cancer_down, commom_colits_down)
  #only cancer
  only_can_up = setdiff(commom_cancer_up,commom_all_up)
  only_can_down = setdiff(commom_cancer_down, commom_all_down)
  #only colits
  only_col_up = setdiff(commom_colits_up,commom_all_up)
  only_col_down = setdiff(commom_colits_down, commom_all_down)
  
  #statistic table
  diff_genes = data.frame(up=1:9, down=1:9, row.names=c("2week", "4week","7week", "10week","Colits_2and4", "Cancer_7and10", "all_week", "only_cancer", "only_colits"))
  diff_genes["2week",] = c(length(diff_genes_2w[[1]]), length(diff_genes_2w[[2]]))
  diff_genes["4week",] = c(length(diff_genes_4w[[1]]), length(diff_genes_4w[[2]]))
  diff_genes["7week",] = c(length(diff_genes_7w[[1]]), length(diff_genes_7w[[2]]))
  diff_genes["10week",] = c(length(diff_genes_10w[[1]]), length(diff_genes_10w[[2]]))
  diff_genes["Colits_2and4",] = c(length(commom_colits_up), length(commom_colits_down))
  diff_genes["Cancer_7and10",] = c(length(commom_cancer_up), length(commom_cancer_down))
  diff_genes["all_week",] = c(length(commom_all_up), length(commom_all_down))
  diff_genes["only_cancer",] = c(length(only_can_up), length(only_can_down))
  diff_genes["only_colits",] = c(length(only_col_up), length(only_col_down))
  summary$total_table = diff_genes
  
  summary$commom_colits_list = list(commom_colits_up, commom_colits_down)
  summary$commom_cancer_list = list(commom_cancer_up, commom_cancer_down)
  summary$commom_col_list = list(only_col_up, only_col_down)
  summary$commom_can_list = list(only_can_up, only_can_down)
  summary$commom_all_list = list(commom_all_up, commom_all_down)
  return(summary)
}

summary = Summary_data("H3K27ac", "enhancer", diff_gene_in_mouse, "RNA")
summary_chip = Summary_data("H3K27ac", "enhancer", diff_gene_for_chip, "Chip")




enrichment = function (list1, list2, mark){
  diff_genes = data.frame(matrix(NA,6,6), row.names=c("cancer_up", "cancer_down","colits_up", "colits_down","all_up", "all_down"))
  colnames(diff_genes) = c("state", "shared_number", "RNA_number", "left_gene_number", "chip_gene_number",  "P_value")
  
  #cancer state
  share_cancer_up = intersect(list1$commom_cancer_list[[1]], list2$commom_cancer_list[[1]])
  pValue = phyper(length(share_cancer_up), length(list1$commom_cancer_list[[1]]), 20000 - length(list1$commom_cancer_list[[1]]), length(list2$commom_cancer_list[[1]]), lower.tail = F)
  diff_genes["cancer_up",] = c(mark, length(share_cancer_up), length(list1$commom_cancer_list[[1]]), 20000 - length(list1$commom_cancer_list[[1]]), length(list2$commom_cancer_list[[1]]), pValue)
  share_cancer_down = intersect(list1$commom_cancer_list[[2]], list2$commom_cancer_list[[2]])
  pValue = phyper(length(share_cancer_down), length(list1$commom_cancer_list[[2]]), 20000 - length(list1$commom_cancer_list[[2]]), length(list2$commom_cancer_list[[2]]), lower.tail = F)
  diff_genes["cancer_down",] = c(mark, length(share_cancer_down), length(list1$commom_cancer_list[[2]]), 20000 - length(list1$commom_cancer_list[[2]]), length(list2$commom_cancer_list[[2]]), pValue)
  
  #colits state
  share_colits_up = intersect(list1$commom_colits_list[[1]], list2$commom_colits_list[[1]])
  pValue = phyper(length(share_colits_up), length(list1$commom_colits_list[[1]]), 20000 - length(list1$commom_colits_list[[1]]), length(list2$commom_colits_list[[1]]), lower.tail = F)
  diff_genes["colits_up",] = c(mark, length(share_colits_up), length(list1$commom_colits_list[[1]]), 20000 - length(list1$commom_colits_list[[1]]), length(list2$commom_colits_list[[1]]), pValue)
  share_colits_down = intersect(list1$commom_colits_list[[2]], list2$commom_colits_list[[2]])
  pValue = phyper(length(share_colits_down), length(list1$commom_colits_list[[2]]), 20000 - length(list1$commom_colits_list[[2]]), length(list2$commom_colits_list[[2]]), lower.tail = F)
  diff_genes["colits_down",] = c(mark, length(share_colits_down), length(list1$commom_colits_list[[2]]), 20000 - length(list1$commom_colits_list[[2]]), length(list2$commom_colits_list[[2]]), pValue)

  #all week
  share_all_up = intersect(list1$commom_all_list[[1]], list2$commom_all_list[[1]])
  pValue = phyper(length(share_all_up), length(list1$commom_all_list[[1]]), 20000 - length(list1$commom_all_list[[1]]), length(list2$commom_all_list[[1]]), lower.tail = F)
  diff_genes["all_up",] = c(mark, length(share_all_up), length(list1$commom_all_list[[1]]), 20000 - length(list1$commom_all_list[[1]]), length(list2$commom_all_list[[1]]), pValue)
  share_all_down = intersect(list1$commom_all_list[[2]], list2$commom_all_list[[2]])
  pValue = phyper(length(share_all_down), length(list1$commom_all_list[[2]]), 20000 - length(list1$commom_all_list[[2]]), length(list2$commom_all_list[[2]]), lower.tail = F)
  diff_genes["all_down",] = c(mark, length(share_all_down), length(list1$commom_all_list[[2]]), 20000 - length(list1$commom_all_list[[2]]), length(list2$commom_all_list[[2]]), pValue)
  return(diff_genes)
}



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




write.table(commom_can_up, "commom_can_up.txt", sep="\n", quote=F, row.names=T, col.names = T, na="NA", eol="\n")
write.table(commom_ccan_up, "commom_ccan_up.txt", sep="\n", quote=F, row.names=T, col.names = T, na="NA", eol="\n")

human_orthologs = mouse2human(commom_ccan_up)
cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
up_symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(human_orthologs), columns=cols, keytype="ENSEMBL")
up_list = unique(up_symbol[, "SYMBOL"])

human_orthologs = mouse2human(commom_ccan_down)
cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
down_symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(human_orthologs), columns=cols, keytype="ENSEMBL")
down_list = unique(down_symbol[, "SYMBOL"]) 

CRCgwasList = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/CRCgwasGene.txt",sep = "\t", header=F)
CRCgwasList = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/CRCmappedGenes.txt",sep = "\t", header=F)

CRCgwasList = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(CRCgwasList$V1), columns=cols, keytype="ENSEMBL")
CRCgwasList = unique(CRCgwasList[, "SYMBOL"])
intersect(CRCgwasList, up_list)



####20190822
























###########################################################
#promoter
df_all = read.delim("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/all_diff_data_H3K4me3.txt",sep = "\t", header=T)
#df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis_0418_non_unique_bam/all_diff_data.txt", sep = "\t", header=T)
mouse_gene_transcript = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_mouse_gene_transcript.txt", sep="\t", header=T)
orthologs_table = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_humanGene_mouseGene.txt", sep="\t", header=T)
peak_gene = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/promoter/promoter.peak.unique.ID.gene_name.bed", sep="\t", header=T)


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
  ind = match(peak_list, peak_gene_table$ID)
  human_genes = peak_gene_table[ind,]
  gene_list = unique(na.omit(human_genes$gene_id))
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


print(c(length(diff_genes_10w[[1]]), length(diff_genes_10w[[2]])))
print(c(length(diff_genes_7w[[1]]), length(diff_genes_7w[[2]])))
print(c(length(diff_genes_4w[[1]]), length(diff_genes_4w[[2]])))
print(c(length(diff_genes_2w[[1]]), length(diff_genes_2w[[2]])))

print(c(length(commom_cancer_up), length(commom_cancer_down)))
print(c(length(commom_colits_up), length(commom_colits_down)))
print(c(length(commom_all_up), length(commom_all_down)))


###is this reasonable?
commom_can_up = setdiff(commom_cancer_up,commom_all_up)
commom_can_down = setdiff(commom_cancer_down, commom_all_down)

commom_col_up = setdiff(commom_colits_up,commom_all_up)
commom_col_down = setdiff(commom_colits_down, commom_all_down)

commom_colits_list_chip = list(commom_colits_up, commom_colits_down)
commom_cancer_list_chip = list(commom_cancer_up, commom_cancer_down)
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


commom_colits_list = list(commom_colits_up, commom_colits_down)
commom_cancer_list = list(commom_cancer_up, commom_cancer_down)
commom_col_list = list(commom_col_up, commom_col_down)
commom_can_list = list(commom_can_up, commom_can_down)
commom_all_list = list(commom_all_up, commom_all_down)

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


commom_colits_list = list(commom_colits_up, commom_colits_down)
commom_cancer_list = list(commom_cancer_up, commom_cancer_down)
commom_col_list = list(commom_col_up, commom_col_down)
commom_can_list = list(commom_can_up, commom_can_down)
commom_all_list = list(commom_all_up, commom_all_down)

#cancer
commom_ccan_up = intersect(commom_cancer_list[[1]], commom_cancer_list_chip[[1]])
pValue = phyper(length(commom_ccan_up), length(commom_cancer_list[[1]]), 20000 - length(commom_cancer_list[[1]]), length(commom_cancer_list_chip[[1]]), lower.tail = F)
print(paste("up", length(commom_ccan_up), length(commom_cancer_list[[1]]), 20000 - length(commom_cancer_list[[1]]), length(commom_cancer_list_chip[[1]]), pValue, sep=","))

commom_ccan_down = intersect(commom_cancer_list[[2]], commom_cancer_list_chip[[2]])
pValue = phyper(length(commom_ccan_down), length(commom_cancer_list[[2]]), 20000 - length(commom_cancer_list[[2]]), length(commom_cancer_list_chip[[2]]), lower.tail = F)
print(paste("down", length(commom_ccan_down), length(commom_cancer_list[[2]]), 20000 - length(commom_cancer_list[[2]]), length(commom_cancer_list_chip[[2]]), pValue, sep=","))

#colites

commom = intersect(commom_colits_list[[1]], commom_colits_list_chip[[1]])
pValue = phyper(length(commom), length(commom_colits_list[[1]]), 20000 - length(commom_colits_list[[1]]), length(commom_colits_list_chip[[1]]), lower.tail = F)
print(paste("up", length(commom), length(commom_colits_list[[1]]), 20000 - length(commom_colits_list[[1]]), length(commom_colits_list_chip[[1]]), pValue, sep=","))

commom = intersect(commom_colits_list[[2]], commom_colits_list_chip[[2]])
pValue = phyper(length(commom), length(commom_colits_list[[2]]), 20000 - length(commom_colits_list[[2]]), length(commom_colits_list_chip[[2]]), lower.tail = F)
print(paste("down", length(commom), length(commom_colits_list[[2]]), 20000 - length(commom_colits_list[[2]]), length(commom_colits_list_chip[[2]]), pValue, sep=","))

#promoter:
commom_cancer_promoter = list(commom_ccan_up, commom_ccan_down)
commom_colits_promoter = list(intersect(commom_colits_list[[1]], commom_colits_list_chip[[1]]), intersect(commom_colits_list[[2]], commom_colits_list_chip[[2]]))


###--------------------------gene list
identy1 = intersect(commom_cancer_promoter[[1]], commom_cancer_enhancer[[1]])
print(length(identy1))
identy2 = intersect(commom_cancer_promoter[[2]], commom_cancer_enhancer[[2]])
print(length(identy2))
identy3 = intersect(commom_colits_promoter[[1]], commom_colits_enhancer[[1]])
print(length(identy3))
identy4 = intersect(commom_colits_promoter[[2]], commom_colits_enhancer[[2]])
print(length(identy4))

#write.table(identy1 , "commom_cancer_up.txt", sep="\n", quote=F, row.names=F, col.names = F, na="NA", eol="\n")




cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
identy5 = intersect(identy1, identy3)
up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(identy5), columns=cols, keytype="ENSEMBL")
identy6 = intersect(identy2, identy4)
up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(identy6), columns=cols, keytype="ENSEMBL")


#GWAS gene
path_gwas = "/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/gwas_gene/gene_list.txt"
CRCgwasList = as.vector(read.delim(path_gwas, sep = "\t", header=F)$V1)


library(org.Hs.eg.db)
#gwas promoter
up_cancer_gwas = intersect(CRCgwasList, mouse2human(commom_cancer_promoter[[1]]))
down_cancer_gwas = intersect(CRCgwasList, mouse2human(commom_cancer_promoter[[2]]))
up_colits_gwas = intersect(CRCgwasList, mouse2human(commom_colits_promoter[[1]]))
down_colits_gwas = intersect(CRCgwasList, mouse2human(commom_colits_promoter[[2]]))


print(length(up_cancer_gwas))
symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(up_cancer_gwas), columns=cols, keytype="ENSEMBL")
print(symbol)
print(length(down_cancer_gwas))
print(length(up_colits_gwas))
symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(up_colits_gwas), columns=cols, keytype="ENSEMBL")
print(symbol)
print(length(down_colits_gwas))


#gwas enhancer
up_cancer_gwas = intersect(CRCgwasList, mouse2human(commom_cancer_enhancer[[1]]))
down_cancer_gwas = intersect(CRCgwasList, mouse2human(commom_cancer_enhancer[[2]]))
up_colits_gwas = intersect(CRCgwasList, mouse2human(commom_colits_enhancer[[1]]))
down_colits_gwas = intersect(CRCgwasList, mouse2human(commom_colits_enhancer[[2]]))


print(length(up_cancer_gwas))
symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(up_cancer_gwas), columns=cols, keytype="ENSEMBL")
print(symbol)
print(length(down_cancer_gwas))
symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(down_cancer_gwas), columns=cols, keytype="ENSEMBL")
print(symbol)
print(length(up_colits_gwas))
symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(up_colits_gwas), columns=cols, keytype="ENSEMBL")
print(symbol)
print(length(down_colits_gwas))


###Cmap-------------------------------------------------
df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis/all_diff_data.txt",sep = "\t", header=T)
col_1 = "log2FoldChange.tumorVSctrl"
col_2 = "padj.tumorVSctrl"
week_df = df_all[,c(col_1, col_2)]
cancer_df = week_df[week_df$padj.tumorVSctrl < 0.01,]
can_df = cancer_df[order(cancer_df$log2FoldChange.tumorVSctrl),]
cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(mouse2human(row.names(head(can_df, 400)))), columns=cols, keytype="ENSEMBL")
up_gene = symbol[,"SYMBOL"]
write.table(up_gene , "cancer_up.txt", sep="\n", quote=F, row.names=F, col.names = F, na="NA", eol="\n")
symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(mouse2human(row.names(tail(can_df, 400)))), columns=cols, keytype="ENSEMBL")
down_gene = symbol[,"SYMBOL"]
write.table(down_gene , "cancer_down.txt", sep="\n", quote=F, row.names=F, col.names = F, na="NA", eol="\n")





diff_gene_in_mouse = function(col_1, col_2, fold_change=2 ,id_type="ENSEMBL"){
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

tummor_vs_control = diff_gene_in_mouse(col_1 = "log2FoldChange.tumorVSctrl", col_2 = "padj.tumorVSctrl")
# divide groups
symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(mouse2human(tummor_vs_control[[1]])), columns=cols, keytype="ENSEMBL")

length(symbol)


phyper(4,10,477,254, lower.tail = FALSE)