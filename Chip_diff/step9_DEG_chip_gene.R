#--------------install package-----------------------------------
install_biomaRt = function (){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("biomaRt", version = "3.8")
}


#--------------analysis-------------------------------------------
library("biomaRt")
library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("DT")

mouse_gene_transcript = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_mouse_gene_transcript.txt", sep="\t", header=T)
orthologs_table = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_humanGene_mouseGene.txt", sep="\t", header=T)
gene_detail = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_mouse_geneid_symbol.txt", sep="\t", header=T)
human_gene_detail = read.delim("/home/zhihl/Project/CRC/RNA_analysis/mart_export_human_geneid_symbol.txt", sep="\t", header=T)
#1. export biomart data

export_biomart_geneid_symbol = function(species, name){
  #species: "Mus musculus", "Homo sapiens"
  #name: mouse, human
  select_dataset=paste(tolower(substring(strsplit(species," ")[[1]][1],1,1)),
                       strsplit(species," ")[[1]][2],"_gene_ensembl",sep = "")
  ensembl=useMart("ensembl", dataset=select_dataset)
  gene_detail=getBM(attributes=c("ensembl_gene_id", "external_gene_name"),mart = ensembl)
  write.table(gene_detail, file = paste("/home/zhihl/Project/CRC/RNA_analysis/mart_export_", name, "_geneid_symbol.txt", sep=""), row.names = T, sep="\t", col.names=T)
}

#export_biomart_geneid_symbol("Homo sapiens", "human")

#2. convert essemble to entrez 
essemble_to_entrez = function (deg_list){
  entrez_map = as.data.frame(unlist(as.list(org.Mm.egENSEMBL2EG)))
  ind<-match(deg_list, rownames(entrez_map))
  deg_eg = entrez_map[ind,]
  return(deg_eg)
}

#3. convert ensembl to symbol
ensembl_to_symbol = function(ensemble_list, name){
  if (name == "mouse")
  {symbol = gene_detail[match(ensemble_list, gene_detail$ensembl_gene_id),"external_gene_name"]}
  if (name == "human")
  {symbol = human_gene_detail[match(ensemble_list, human_gene_detail$ensembl_gene_id),"external_gene_name"]}
  return(as.vector(symbol)) 
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

human2mouse = function(deg_list){
  orthologs = orthologs_table
  ind = match(deg_list, orthologs$Gene.stable.ID)
  human_genes = orthologs[ind,]
  gene_list = unique(na.omit(human_genes[,3]))
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



#6. for Chip
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
  
  if (id_type == "ENSEMBL"){
    return(list(unique(up_gene), unique(down_gene)))
  }
  
  cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
  up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(up_gene), columns=cols, keytype="ENSEMBL")
  up_list = unique(up_symbol[, id_type])
  down_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(down_gene), columns=cols, keytype="ENSEMBL")
  down_list = unique(down_symbol[, id_type])
  return (list(up_list, down_list))
}

#7. for RNA
diff_gene_in_mouse = function(df_all, col_1, col_2, fold_change=1 ,id_type="ENSEMBL"){
  #col_1: col name of log2 fold change
  #col_2: col name of padj
  #fold_change: log2 fold change
  #id_type: "SYMBOL", "ENTRIZEID" ....
  #return list of up and down genes
  
  week_df = df_all[,c(col_1, col_2)]
  diff_week_up = rownames(na.omit(week_df[week_df[,1] > fold_change & week_df[,2] < 0.01, ]))
  diff_week_down = rownames(na.omit(week_df[week_df[,1] < -(fold_change) & week_df[,2] < 0.01, ]))
  
  
  if (id_type == "ENSEMBL"){
    return(list(unique(diff_week_up), unique(diff_week_down)))
  }
  
  cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
  up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(diff_week_up), columns=cols, keytype="ENSEMBL")
  up_list = unique(up_symbol[, id_type])
  down_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(diff_week_down), columns=cols, keytype="ENSEMBL")
  down_list = unique(down_symbol[, id_type])
  return (list(up_list, down_list))
}

#8. statistic
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
  
  
  #tummor_vs_control = diff_gene_in_mouse(df_rna, "log2FoldChange.tumorVSctrl", "padj.tumorVSctrl", 1, "ENSEMBL")
  
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

#9.enrichment RNA and Chip
enrichment = function (list1, list2, mark){
  diff_genes = data.frame(matrix(NA,6,6), row.names=c("cancer_up", "cancer_down","colits_up", "colits_down","all_up", "all_down"))
  colnames(diff_genes) = c("marker", "shared_number", "RNA_number", "left_gene_number", "chip_gene_number",  "P_value")
  
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



#10. two mark result
two_mark = function(summary_enhancer, summary_promoter, summary_rna, marks){
  df = data.frame(matrix(NA,18,6))
  
  #two_mark, state, up_down, geneID, symbol, withRNA
  
  ####enhancer, promoter
  
  #cancer 
  up_enhancer_promoter = intersect(summary_enhancer$commom_cancer_list[[1]], summary_promoter$commom_cancer_list[[1]])
  withRNA = intersect(up_enhancer_promoter, summary_rna$commom_cancer_list[[1]])
  
  df_1 =  data.frame(two_mark = rep(marks, length(up_enhancer_promoter)), state = rep("cancer", length(up_enhancer_promoter)),
  up_down = rep("up", length(up_enhancer_promoter)), geneID = up_enhancer_promoter, 
  symbol = gene_detail[match(up_enhancer_promoter, gene_detail$ensembl_gene_id),"external_gene_name"], withRNA = ifelse(up_enhancer_promoter%in%withRNA,"true","false"))
  
  down_enhancer_promoter = intersect(summary_enhancer$commom_cancer_list[[2]], summary_promoter$commom_cancer_list[[2]])
  withRNA = intersect(down_enhancer_promoter, summary_rna$commom_cancer_list[[2]])
  df_2 =  data.frame(two_mark = rep(marks, length(down_enhancer_promoter)), state = rep("cancer", length(down_enhancer_promoter)),
                     up_down = rep("down", length(down_enhancer_promoter)), geneID = down_enhancer_promoter, 
                     symbol = gene_detail[match(down_enhancer_promoter, gene_detail$ensembl_gene_id),"external_gene_name"], withRNA = ifelse(down_enhancer_promoter%in%withRNA,"true","false"))
  new_df = rbind(df_1, df_2)
  #colits
  up_enhancer_promoter = intersect(summary_enhancer$commom_colits_list[[1]], summary_promoter$commom_colits_list[[1]])
  withRNA = intersect(up_enhancer_promoter, summary_rna$commom_colits_list[[1]])
  
  df_3 =  data.frame(two_mark = rep(marks, length(up_enhancer_promoter)), state = rep("colits", length(up_enhancer_promoter)),
                     up_down = rep("up", length(up_enhancer_promoter)), geneID = up_enhancer_promoter, 
                     symbol = gene_detail[match(up_enhancer_promoter, gene_detail$ensembl_gene_id),"external_gene_name"], withRNA = ifelse(up_enhancer_promoter%in%withRNA,"true","false"))
  
  new_df = rbind(new_df, df_3)
  down_enhancer_promoter = intersect(summary_enhancer$commom_colits_list[[2]], summary_promoter$commom_colits_list[[2]])
  withRNA = intersect(down_enhancer_promoter, summary_rna$commom_colits_list[[2]])
  df_4 =  data.frame(two_mark = rep(marks, length(down_enhancer_promoter)), state = rep("colits", length(down_enhancer_promoter)),
                     up_down = rep("down", length(down_enhancer_promoter)), geneID = down_enhancer_promoter, 
                     symbol = gene_detail[match(down_enhancer_promoter, gene_detail$ensembl_gene_id),"external_gene_name"], withRNA = ifelse(down_enhancer_promoter%in%withRNA,"true","false"))
  new_df = rbind(new_df, df_4)
  
  #all
  up_enhancer_promoter = intersect(summary_enhancer$commom_all_list[[1]], summary_promoter$commom_all_list[[1]])
  withRNA = intersect(up_enhancer_promoter, summary_rna$commom_all_list[[1]])
  
  df_5 =  data.frame(two_mark = rep(marks, length(up_enhancer_promoter)), state = rep("all", length(up_enhancer_promoter)),
                     up_down = rep("up", length(up_enhancer_promoter)), geneID = up_enhancer_promoter, 
                     symbol = gene_detail[match(up_enhancer_promoter, gene_detail$ensembl_gene_id),"external_gene_name"], withRNA = ifelse(up_enhancer_promoter%in%withRNA,"true","false"))
  new_df = rbind(new_df, df_5)
  down_enhancer_promoter = intersect(summary_enhancer$commom_all_list[[2]], summary_promoter$commom_all_list[[2]])
  withRNA = intersect(down_enhancer_promoter, summary_rna$commom_all_list[[2]])
  df_6 =  data.frame(two_mark = rep(marks, length(down_enhancer_promoter)), state = rep("all", length(down_enhancer_promoter)),
                     up_down = rep("down", length(down_enhancer_promoter)), geneID = down_enhancer_promoter, 
                     symbol = gene_detail[match(down_enhancer_promoter, gene_detail$ensembl_gene_id),"external_gene_name"], withRNA = ifelse(down_enhancer_promoter%in%withRNA,"true","false"))
  new_df = rbind(new_df, df_6)
  
  
  
  
  
  
  
  
  return(new_df)
}

#11. reverse list
reverse_list = function(summary_list){
  summary_reversed = summary_list
  summary_reversed$commom_cancer_list[[2]] = summary_list$commom_cancer_list[[1]]
  summary_reversed$commom_cancer_list[[1]] = summary_list$commom_cancer_list[[2]]
  summary_reversed$commom_colits_list[[2]] = summary_list$commom_colits_list[[1]]
  summary_reversed$commom_colits_list[[1]] = summary_list$commom_cancer_list[[2]]
  summary_reversed$commom_all_list[[2]] = summary_list$commom_all_list[[1]]
  summary_reversed$commom_all_list[[1]] = summary_list$commom_all_list[[2]]
  summary_reversed$commom_can_list[[2]] = summary_list$commom_can_list[[1]]
  summary_reversed$commom_can_list[[1]] = summary_list$commom_can_list[[2]]
  summary_reversed$commom_col_list[[2]] = summary_list$commom_col_list[[1]]
  summary_reversed$commom_col_list[[1]] = summary_list$commom_col_list[[2]]
  return(summary_reversed)
}

#12. gwas gene
gwas_gene = function(summary_rna, name){
  #summary_rna, summary_enhancer, summary_promoter, summary_repressed, summary_heterochromatin
  #name: "rna", "enhancer", "promoter", "repressed", "heterochromatin"
  path_gwas = "/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/gwas_gene/gene_list.txt"
  CRCgwasList = as.vector(read.delim(path_gwas, sep = "\t", header=F)$V1)
  
  #mark, state, up_down, mouseGeneID, mouseSymbol, humanGeneID, humanSymbol
  #cancer
  up_cancer_gwas = intersect(CRCgwasList, mouse2human(summary_rna$commom_cancer_list[[1]]))
  length_num = length(up_cancer_gwas)
  df_1 = data.frame(GWAS=rep("GWAS", length_num), mark=rep(name, length_num), 
                    state=rep("cancer", length_num), up_down=rep("up", length_num), 
                    mouseGeneID=human2mouse(up_cancer_gwas),
                    mouseSymbol=ensembl_to_symbol(human2mouse(up_cancer_gwas), "mouse"), 
                    humanGeneID=up_cancer_gwas,  humanSymbol=ensembl_to_symbol(up_cancer_gwas, "human"))
  
  down_cancer_gwas = intersect(CRCgwasList, mouse2human(summary_rna$commom_cancer_list[[2]]))
  length_num = length(down_cancer_gwas)
  df_2 = data.frame(GWAS=rep("GWAS", length_num), mark=rep(name, length_num), 
                    state=rep("cancer", length_num), up_down=rep("down", length_num),
                    mouseGeneID=human2mouse(down_cancer_gwas),
                    mouseSymbol=ensembl_to_symbol(human2mouse(down_cancer_gwas), "mouse"), 
                    humanGeneID=down_cancer_gwas, humanSymbol=ensembl_to_symbol(down_cancer_gwas, "human"))
  
  new_df = rbind(df_1, df_2)
  #colits
  #up_colits_gwas = intersect(CRCgwasList, mouse2human(summary_rna$commom_colits_list[[1]]))
  up_colits_gwas = intersect(CRCgwasList, mouse2human(summary_rna$commom_colits_list[[1]]))
  length_num = length(up_colits_gwas)
  df_3 = data.frame(GWAS=rep("GWAS", length_num), mark=rep(name, length_num), 
                    state=rep("colits", length_num), up_down=rep("up", length_num),
                    mouseGeneID=human2mouse(up_colits_gwas),
                    mouseSymbol=ensembl_to_symbol(human2mouse(up_colits_gwas), "mouse"), 
                    humanGeneID=up_colits_gwas, humanSymbol=ensembl_to_symbol(up_colits_gwas, "human"))
  new_df = rbind(new_df, df_3)
  
  down_colits_gwas = intersect(CRCgwasList, mouse2human(summary_rna$commom_colits_list[[2]]))
  length_num = length(down_colits_gwas)
  df_4 = data.frame(GWAS=rep("GWAS", length_num), mark=rep(name, length_num), 
                    state=rep("colits", length_num), up_down=rep("down", length_num),
                    mouseGeneID=human2mouse(down_colits_gwas),
                    mouseSymbol=ensembl_to_symbol(human2mouse(down_colits_gwas), "mouse"), 
                    humanGeneID=down_colits_gwas, humanSymbol=ensembl_to_symbol(down_colits_gwas, "human"))
  new_df = rbind(new_df, df_4)
  
  #all
  up_all_gwas = intersect(CRCgwasList, mouse2human(summary_rna$commom_all_list[[1]]))
  length_num = length(up_all_gwas)
  df_5 = data.frame(GWAS=rep("GWAS", length_num), mark=rep(name, length_num), 
                    state=rep("all", length_num), up_down=rep("up", length_num),
                    mouseGeneID=human2mouse(up_all_gwas),
                    mouseSymbol=ensembl_to_symbol(human2mouse(up_all_gwas), "mouse"), 
                    humanGeneID=up_all_gwas, humanSymbol=ensembl_to_symbol(up_all_gwas, "human"))
  new_df = rbind(new_df, df_5)
  
  down_all_gwas = intersect(CRCgwasList, mouse2human(summary_rna$commom_all_list[[2]]))
  length_num = length(down_all_gwas)
  df_6 = data.frame(GWAS=rep("GWAS", length_num), mark=rep(name, length_num), 
                    state=rep("all", length_num), up_down=rep("down", length_num),
                    mouseGeneID=human2mouse(down_all_gwas),
                    mouseSymbol=ensembl_to_symbol(human2mouse(down_all_gwas), "mouse"), 
                    humanGeneID=down_all_gwas, humanSymbol=ensembl_to_symbol(down_all_gwas, "human"))
  
  new_df = rbind(new_df, df_6)
  return(new_df)
}

#13. add gene symbol to source data
add_symbol = function(table_path, type, out_path, peak_gene_path){
  df_all = read.delim(table_path,sep = "\t", header=T)
  if (type == "RNA"){
    ensemblID = rownames(df_all)
    symbol = ensembl_to_symbol(ensemblID, "mouse")
    new_df = cbind.data.frame(ensemblID, df_all)
    new_df = cbind.data.frame(symbol, new_df)
  }
  if (type == "Chip"){
    peak_gene = read.delim(peak_gene_path, sep="\t", header=T)
    peak_gene_table = peak_gene
    peak_id = rownames(df_all)
    ind = match(peak_id, peak_gene_table$ID)
    ensemblID = peak_gene_table[ind, "gene_id"]
    symbol = ensembl_to_symbol(ensemblID, "mouse")
    new_df = cbind.data.frame(peak_id, df_all)
    new_df = cbind.data.frame(ensemblID, new_df)
    new_df = cbind.data.frame(symbol, new_df)
  }
  
  write.table(new_df, out_path, sep="\t", quote=F, row.names=F, col.names = T, na="NA", eol="\n")
}
add_symbol("/home/zhihl/Project/CRC/RNA_analysis/all_diff_data.txt", "RNA", "/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/final_result/all_diff_data_add_symbol.txt", "")
add_symbol("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/all_diff_data_H3K27ac.txt", "Chip",
           "/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/final_result/all_diff_data_H3K27ac_add_symbol.txt",
           "/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/enhancer/enhancer.peak.unique.ID.gene_name.bed")
add_symbol("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/all_diff_data_H3K4me3.txt", "Chip",
           "/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/final_result/all_diff_data_H3K4me3_add_symbol.txt",
           "/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/promoter/promoter.peak.unique.ID.gene_name.bed")
add_symbol("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/all_diff_data_H3K27me3.txt", "Chip",
           "/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/final_result/all_diff_data_H3K27me3_add_symbol.txt",
           "/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/repressed/repressed.peak.unique.ID.gene_name.bed")
add_symbol("/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/all_diff_data_H3K9me3.txt", "Chip",
           "/home/zhihl/Project/CRC/Chip_analysis/dff/version0821/final_result/all_diff_data_H3K9me3_add_symbol.txt",
           "/home/zhihl/Project/CRC/Chip_analysis/peak_dir_0821/heterochromatin/heterochromatin.peak.unique.ID.gene_name.bed")



#-------------------------other GWAS data-------------------------------------------------------------------------------------
#human_orthologs = mouse2human(commom_ccan_up)
#cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
#up_symbol = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(human_orthologs), columns=cols, keytype="ENSEMBL")
#up_list = unique(up_symbol[, "SYMBOL"])


#CRCgwasList = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/CRCgwasGene.txt",sep = "\t", header=F)
#CRCgwasList = read.delim("/home/zhihl/Project/CRC/Chip_analysis/peak_gene_mapping/CRCmappedGenes.txt",sep = "\t", header=F)

#CRCgwasList = AnnotationDbi::select(org.Hs.eg.db, keys=as.vector(CRCgwasList$V1), columns=cols, keytype="ENSEMBL")
#CRCgwasList = unique(CRCgwasList[, "SYMBOL"])
#intersect(CRCgwasList, up_list)


#--------------------------translate-------------------------------------------------------------------------------------------
# translate to symbol
#cols <- c("SYMBOL", "ENTREZID", "ENSEMBL")
#up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(identy5), columns=cols, keytype="ENSEMBL")
#up_symbol = AnnotationDbi::select(org.Mm.eg.db, keys=as.vector(identy6), columns=cols, keytype="ENSEMBL")
