#--------------install package-----------------------------------
install_biomaRt = function (){
  if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  BiocManager::install("biomaRt", version = "3.8")
}

#--------------analysis-------------------------------------------
library("biomaRt")
library("org.Hs.eg.db")

df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis/all_diff_data.txt",sep = "\t", header=T)
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
diff_gene = function(col_1, col_2, fold_change=2 ,id_type="SYMBOL"){
  #col_1: col name of log2 fold change
  #col_2: col name of padj
  #fold_change: log2 fold change
  #id_type: "SYMBOL", "ENTRIZEID" ....
  #return list of up and down genes
  
  week_df = df_all[,c(col_1, col_2)]
  diff_week_up = rownames(na.omit(week_df[week_df[,1] > fold_change & week_df[,2] < 0.01, ]))
  diff_week_down = rownames(na.omit(week_df[week_df[,1] < -(fold_change) & week_df[,2] < 0.01, ]))
  down_gene = mouse2human(diff_week_down)
  up_gene = mouse2human(diff_week_up)
  
  cols <- c("SYMBOL", "ENTREZID")
  up_symbol = select(org.Hs.eg.db, keys=as.vector(up_gene), columns=cols, keytype="ENSEMBL")
  up_list = unique(up_symbol[, id_type])
  down_symbol = select(org.Hs.eg.db, keys=as.vector(down_gene), columns=cols, keytype="ENSEMBL")
  down_list = unique(down_symbol[, id_type])
  return (list(up_list, down_list))
}

diff_genes_10w = diff_gene("log2FoldChange.10wkVSctrl", "padj.10wkVSctrl")
diff_genes_7w = diff_gene("log2FoldChange.7wkVSctrl", "padj.7wkVSctrl")
diff_genes_4w = diff_gene("log2FoldChange.4wkVSctrl", "padj.4wkVSctrl")
diff_genes_2w = diff_gene("log2FoldChange.2wkVSctrl", "padj.2wkVSctrl")


#paper 2
up_genes = read.delim("/home/zhihl/Project/CRC/RNA_analysis/reference_diff_gene_paper1/paper2_up.txt",sep = "\t", header=F)
up_gene_unique = unique(as.vector(up_genes[,1]))
commom = intersect(diff_genes_10w[[1]], up_gene_unique)


down_genes = read.delim("/home/zhihl/Project/CRC/RNA_analysis/reference_diff_gene_paper1/paper2_down.txt",sep = "\t", header=F)
down_gene_unique = unique(as.vector(down_genes[,1]))
commom = intersect(diff_genes_10w[[2]], down_gene_unique)


#paper 3:Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value
f_intersect = function (x, diff_genes_w){
  total_gene = 20000
  data_1 = data_all[data_all[, x[1]] > 0 & data_all[, x[2]] < 0.05,]
  up_list = unique(as.vector(data_1$Gene.Symbol))
  #print(up_list)
  data_1 = data_all[data_all[, x[1]] < 0 & data_all[,x[2]] < 0.05,]
  down_list = unique(as.vector(data_1$Gene.Symbol))
  
  print("up gene overlap:")
  commom = intersect(diff_genes_w[[1]], up_list)
  #print(commom)
  #hypergenomic test
  pValue = phyper(length(commom), length(up_list), 20000 - length(up_list), length(diff_genes_w[[1]]), lower.tail = F)
  print(paste(x[1], "up", length(commom), length(up_list), 20000 - length(up_list), length(diff_genes_w[[1]]), pValue, sep=","))
  
  
  print("down gene overlap")
  commom = intersect(diff_genes_w[[2]], down_list)
  #print(commom)
  #hypergenomic test
  pValue = phyper(length(commom), length(down_list), 20000 - length(down_list), length(diff_genes_w[[2]]), lower.tail = F)
  print(paste(x[2], "down", length(commom), length(down_list), 20000 - length(down_list), length(diff_genes_w[[2]]), pValue, sep=","))
  
}

paper3_analysis = function(diff_genes_w){
  data_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis/reference_diff_gene_paper1/paper3/data.txt",sep = "\t", header=T)
  select_subtype_columns = c("logFC.C1vsCx", "logFC.C2vsCx", "logFC.C3vsCx", "logFC.C4vsCx", "logFC.C5vsCx", "logFC.C6vsCx")
  select_adj_columns = c("adjpv.C1vsCx", "adjpv.C2vsCx", "adjpv.C3vsCx", "adjpv.C4vsCx", "adjpv.C5vsCx", "adjpv.C6vsCx")
  select_frame = data.frame(logFC = select_subtype_columns, adjPv = select_adj_columns)
  apply(select_frame, 1, f_intersect, diff_genes_w)
}

result = paper3_analysis(diff_genes_w = diff_genes_10w)
result = paper3_analysis(diff_genes_w = diff_genes_7w)
result = paper3_analysis(diff_genes_w = diff_genes_4w)
result = paper3_analysis(diff_genes_w = diff_genes_2w)


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


result = paper3_analysis(diff_genes_w = commom_all_list)
result = paper3_analysis(diff_genes_w = commom_can_list)
result = paper3_analysis(diff_genes_w = commom_col_list)


#C2
data_1 = data_all[data_all$logFC.C2vsCx > 0 & data_all$adjpv.C2vsCx <0.05,]
up_list = unique(as.vector(data_1$Gene.Symbol))
data_1 = data_all[data_all$logFC.C2vsCx < 0 & data_all$adjpv.C2vsCx <0.05,]
down_list = unique(as.vector(data_1$Gene.Symbol))

commom = intersect(diff_genes_10w[[1]], up_list)
commom
commom = intersect(diff_genes_10w[[2]], down_list)
commom
