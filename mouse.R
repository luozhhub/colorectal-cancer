library("org.Mm.eg.db")

cols <- c("SYMBOL", "ENTREZID")
df_all = read.delim("/home/zhihl/Project/CRC/RNA_analysis/all_diff_data.txt",sep = "\t", header=T)

#add gene symbol in dataframe
add_symbol = function(){
gene_symbol = select(org.Mm.eg.db, keys=as.vector(rownames(df_all)), columns=cols, keytype="ENSEMBL")
gene_symbol_new = na.omit(gene_symbol)
df_add = merge(x = df_all, y = gene_symbol, by.x = "row.names", by.y="ENSEMBL" , all = TRUE)
write.table(df_add, "/home/zhihl/Project/CRC/RNA_analysis/all_diff_data_addSymbol.txt", sep="\t", quote=F, row.names=T, col.names = T, na="NA", eol="\n")
}


#select data
diff_gene_in_mouse = function(col_1, col_2, fold_change=1 ,id_type="SYMBOL"){
  #col_1: col name of log2 fold change
  #col_2: col name of padj
  #fold_change: log2 fold change
  #id_type: "SYMBOL", "ENTRIZEID" ....
  #return list of up and down genes
  
  week_df = df_all[,c(col_1, col_2)]
  diff_week_up = rownames(na.omit(week_df[week_df[,1] > fold_change & week_df[,2] < 0.01, ]))
  diff_week_down = rownames(na.omit(week_df[week_df[,1] < -(fold_change) & week_df[,2] < 0.01, ]))
  
  
  cols <- c("SYMBOL", "ENTREZID")
  up_symbol = select(org.Mm.eg.db, keys=as.vector(diff_week_up), columns=cols, keytype="ENSEMBL")
  up_list = unique(up_symbol[, id_type])
  down_symbol = select(org.Mm.eg.db, keys=as.vector(diff_week_down), columns=cols, keytype="ENSEMBL")
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



#paper 7:
paper7_analysis = function(list_of_genes ){
  p_data = read.delim("/home/zhihl/Project/CRC/RNA_analysis/reference_diff_gene_paper1/paper7/up.txt", sep = "\t", header=T)
  whole_genes = unique(list_of_genes[[1]])
  genes_type1 = unique(p_data$tracking_id)
  commom = intersect(whole_genes, genes_type1)
  pValue = phyper(length(commom), length(genes_type1), 20000 - length(genes_type1), length(whole_genes), lower.tail = F)
  print(paste("up", length(commom), length(genes_type1), 20000 - length(genes_type1), length(whole_genes), pValue, sep=","))
  
  
  p_data = read.delim("/home/zhihl/Project/CRC/RNA_analysis/reference_diff_gene_paper1/paper7/down.txt", sep = "\t", header=T)
  whole_genes = unique(list_of_genes[[2]])
  genes_type1 = unique(p_data$tracking_id)
  commom = intersect(whole_genes, genes_type1)
  pValue = phyper(length(commom), length(genes_type1), 20000 - length(genes_type1), length(whole_genes), lower.tail = F)
  print(paste("down", length(commom), length(genes_type1), 20000 - length(genes_type1), length(whole_genes), pValue, sep=","))
  
}

paper7_analysis(list_of_genes = commom_all_list)
paper7_analysis(list_of_genes = commom_can_list)
paper7_analysis(list_of_genes = commom_col_list)


#paper 8:Neufert C, Becker C, Türeci Ö, et al. Tumor fibroblast–derived epiregulin promotes growth of colitis-associated neoplasms through ERK[J]. The Journal of clinical investigation, 2013, 123(4): 1428-1443.
paper8_analysis = function(list_of_genes ){
  p_data = read.delim("/home/zhihl/Project/CRC/RNA_analysis/reference_diff_gene_paper1/mouse_paper1/data.txt", sep = "\t", header=T)
  whole_genes = unique(list_of_genes[[1]])
  #genes_type1 = unique(p_data$Gene.Symbol[p_data$correctedPvalue..sporCRC..Vs..sporCRC.contr.. < 0.01 & p_data$regulation..sporCRC..Vs..sporCRC.contr.. == "up"])
  genes_type1 = unique(p_data$Gene.Symbol[p_data$correctedPvalue..CAC..Vs..CAC.contr.. < 0.01 & p_data$regulation..CAC..Vs..CAC.contr.. == "up"])
  
  commom = intersect(whole_genes, genes_type1)
  pValue = phyper(length(commom), length(genes_type1), 20000 - length(genes_type1), length(whole_genes), lower.tail = F)
  print(paste("up", length(commom), length(genes_type1), 20000 - length(genes_type1), length(whole_genes), pValue, sep=","))
  
  
  p_data = read.delim("/home/zhihl/Project/CRC/RNA_analysis/reference_diff_gene_paper1/mouse_paper1/data.txt", sep = "\t", header=T)
  whole_genes = unique(list_of_genes[[2]])
  #genes_type1 = unique(p_data$Gene.Symbol[p_data$correctedPvalue..sporCRC..Vs..sporCRC.contr.. < 0.01 & p_data$regulation..sporCRC..Vs..sporCRC.contr.. == "down"])
  genes_type1 = unique(p_data$Gene.Symbol[p_data$correctedPvalue..CAC..Vs..CAC.contr.. < 0.01 & p_data$regulation..CAC..Vs..CAC.contr.. == "down"])
  commom = intersect(whole_genes, genes_type1)
  pValue = phyper(length(commom), length(genes_type1), 20000 - length(genes_type1), length(whole_genes), lower.tail = F)
  print(paste("down", length(commom), length(genes_type1), 20000 - length(genes_type1), length(whole_genes), pValue, sep=","))
  
}

paper8_analysis(list_of_genes = commom_all_list)
paper8_analysis(list_of_genes = commom_can_list)
paper8_analysis(list_of_genes = commom_col_list)


paper8_analysis(list_of_genes = diff_genes_10w)
paper8_analysis(list_of_genes = diff_genes_4w)

paper8_analysis(list_of_genes = tummor_vs_control)

