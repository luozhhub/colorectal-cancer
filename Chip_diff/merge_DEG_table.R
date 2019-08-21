#!/usr/bin/R
args=commandArgs(T)
read_count_table_path = args[1]
output_dir = args[2]
output_name = args[3]
library(DESeq2)

#1. read count table file,return a deseq calss
preprocess = function (count_table_path, replicate_number=3){
  d.raw <- read.delim(count_table_path, sep = "\t")
  d <- d.raw[rowSums(d.raw>3) > 2,]
  grp <- c(rep("control",replicate_number),rep("2Week",replicate_number),rep("4Week",replicate_number),rep("7Week",replicate_number),rep("10Week",replicate_number))
  cData <- data.frame(Time = factor(grp, levels = c("control", "2Week", "4Week", "7Week", "10Week")))
  rownames(cData) <- colnames(d)
  d.deseq <- DESeqDataSetFromMatrix(countData = d, colData = cData,design = ~Time)
  d.deseq <- DESeq(d.deseq)
  return(d.deseq)
}

#1.1 group as "control", "inflam", "tumor"
preprocess_2 = function (count_table_path, replicate_number=3 ){
  d.raw <- read.delim(count_table_path, sep = "\t")
  d <- d.raw[rowSums(d.raw>3) > 2,]
  grp <- c(rep("control",replicate_number),rep("inflam",replicate_number),rep("inflam",replicate_number),rep("tumor",replicate_number),rep("tumor",replicate_number))
  cData <- data.frame(Time = factor(grp, levels = c("control", "inflam", "tumor")))
  rownames(cData) <- colnames(d)
  d.deseq <- DESeqDataSetFromMatrix(countData = d, colData = cData,design = ~Time)
  d.deseq <- DESeq(d.deseq)
  return(d.deseq)
}

#construct the total data
select_data = function(count_table_path){
  data = preprocess(count_table_path = count_table_path, replicate_number = 2)
  data2 = preprocess_2(count_table_path = count_table_path, replicate_number = 2)
  
  ###data select
  group = "Time_2Week_vs_control"
  res_dds_main=na.omit(as.data.frame(results(data, name=group, cooksCutoff=FALSE)))
  deg_2week_VS_ctrl = res_dds_main[,c("log2FoldChange", "padj")]

  group = "Time_4Week_vs_control"
  res_dds_main=na.omit(as.data.frame(results(data, name=group, cooksCutoff=FALSE)))
  deg_4week_VS_ctrl = res_dds_main[,c("log2FoldChange", "padj")]

  group = "Time_7Week_vs_control"
  res_dds_main=na.omit(as.data.frame(results(data, name=group, cooksCutoff=FALSE)))
  deg_7week_VS_ctrl = res_dds_main[,c("log2FoldChange", "padj")]

  group = "Time_10Week_vs_control"
  res_dds_main=na.omit(as.data.frame(results(data, name=group, cooksCutoff=FALSE)))
  deg_10week_VS_ctrl = res_dds_main[,c("log2FoldChange", "padj")]

  contrast = c("Time", "2Week", "4Week")
  res_dds_main=na.omit(as.data.frame(results(data, contrast = contrast, cooksCutoff=FALSE)))
  deg_2week_VS_4week = res_dds_main[,c("log2FoldChange", "padj")]
  contrast = c("Time", "4Week", "7Week")
  res_dds_main=na.omit(as.data.frame(results(data, contrast = contrast, cooksCutoff=FALSE)))
  deg_4week_VS_7week = res_dds_main[,c("log2FoldChange", "padj")]
  contrast = c("Time", "7Week", "10Week")
  res_dds_main=na.omit(as.data.frame(results(data, contrast = contrast, cooksCutoff=FALSE)))
  deg_7week_VS_10week = res_dds_main[,c("log2FoldChange", "padj")]


  ###data2 selecting
  contrast = c("Time", "inflam", "control")
  res_dds_main=na.omit(as.data.frame(results(data2, contrast = contrast, cooksCutoff=FALSE)))
  deg_inflam_VS_ctrl = res_dds_main[,c("log2FoldChange", "padj")]
  contrast = c("Time", "tumor", "control")
  res_dds_main=na.omit(as.data.frame(results(data2, contrast = contrast, cooksCutoff=FALSE)))
  deg_tumor_VS_ctrl = res_dds_main[,c("log2FoldChange", "padj")]
  contrast = c("Time", "tumor", "inflam")
  res_dds_main=na.omit(as.data.frame(results(data2, contrast = contrast, cooksCutoff=FALSE)))
  deg_tumor_VS_inflam = res_dds_main[,c("log2FoldChange", "padj")]


  ###merge data
  #merge 4
  merge_data = merge(x = deg_2week_VS_ctrl, y = deg_4week_VS_ctrl, by = "row.names", all = TRUE, suffixes = c(".2wkVSctrl",".4wkVSctrl"))
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  #merge 7
  merge_data = merge(x = merge_data, y = deg_7week_VS_ctrl, by = "row.names", all = TRUE)
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  #merge 10
  merge_data = merge(x = merge_data, y = deg_10week_VS_ctrl, by = "row.names", all = TRUE, suffixes = c(".7wkVSctrl",".10wkVSctrl"))
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  #merge 2_4
  merge_data = merge(x = merge_data, y = deg_2week_VS_4week, by = "row.names", all = TRUE)
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  
  merge_data = merge(x = merge_data, y = deg_4week_VS_7week, by = "row.names", all = TRUE, suffixes = c(".2wkVS4wk",".4wkVS7wk"))
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  
  merge_data = merge(x = merge_data, y = deg_7week_VS_10week, by = "row.names", all = TRUE)
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  
  merge_data = merge(x = merge_data, y = deg_inflam_VS_ctrl, by = "row.names", all = TRUE, suffixes = c(".7wkVS10wk",".inflVSctrl"))
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  
  merge_data = merge(x = merge_data, y = deg_tumor_VS_ctrl, by = "row.names", all = TRUE)
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  
  merge_data = merge(x = merge_data, y = deg_tumor_VS_inflam, by = "row.names", all = TRUE, suffixes = c(".tumorVSctrl", ".tumorVSinfl"))
  rownames(merge_data) = merge_data$Row.names
  merge_data <- subset(merge_data, select = -Row.names)
  

  return(merge_data)
}

#merge DEG columns, write table
create_total_table = function(read_count_table_path, output_name,  output_dir){
  
  merge_data = select_data(count_table_path = read_count_table_path)
  #merge_data = select_data(count_table_path = "/home/zhihl/Project/CRC/RNA_analysis/count_table_average.txt")
  if (!file_test("-d", output_dir)){ dir.create(output_dir)}
  out_file_path = paste(output_dir, output_name, sep = "/")
  write.table(merge_data, out_file_path, sep="\t", quote=F, row.names=T, col.names = T, na="NA", eol="\n")
}

create_total_table(read_count_table_path = read_count_table_path, output_dir = output_dir, output_name=output_name)
