#set work dir
setwd("/home/zhihl/Project/CRC/RNA_analysis")
#include package, thses packages must be installed before
library(DESeq2)
library("org.Mm.eg.db")
library(ggplot2)


#1. process the read count data, and run deseq2
preprocess = function (){
  d.raw <- read.delim("/home/zhihl/Project/CRC/RNA_analysis/count_table.txt",sep = "\t")
  d <- d.raw[rowSums(d.raw>3) > 2,]
  grp <- c(rep("control",3),rep("2Week",3),rep("4Week",3),rep("7Week",3),rep("10Week",3))
  cData <- data.frame(Time = factor(grp, levels = c("control", "2Week", "4Week", "7Week", "10Week")))
  rownames(cData) <- colnames(d)
  d.deseq <- DESeqDataSetFromMatrix(countData = d, colData = cData,design = ~Time)
  d.deseq <- DESeq(d.deseq)
  return(d.deseq)
}

preprocess_2 = function (){
  d.raw <- read.delim("count_table.txt",sep = "\t")
  d <- d.raw[rowSums(d.raw>3) > 2,]
  grp <- c(rep("control",3),rep("inflam",3),rep("inflam",3),rep("tumor",3),rep("tumor",3))
  cData <- data.frame(Time = factor(grp, levels = c("control", "inflam", "tumor")))
  rownames(cData) <- colnames(d)
  d.deseq <- DESeqDataSetFromMatrix(countData = d, colData = cData,design = ~Time)
  d.deseq <- DESeq(d.deseq)
  return(d.deseq)
}

#normalize
#b=counts(data,normalized=T)

#2. parse the different expression gene
#resultsNames: "Intercept","Time_2Week_vs_control", "Time_4Week_vs_control", "Time_7Week_vs_control", "Time_10Week_vs_control"
#contrast eg. : contrast = c("Time", "4Week", "control")
degfunc = function (group="", data, contrast=""){
  if (group != ""){
    res_dds_main=na.omit(as.data.frame(results(data, name=group, cooksCutoff=FALSE)))
    print("111111")
  }
  if (contrast != ""){
    res_dds_main=na.omit(as.data.frame(results(data, contrast = contrast, cooksCutoff=FALSE)))
  }
  degs = res_dds_main[res_dds_main$log2FoldChange > 1 & res_dds_main$padj < 0.05,]
  deg_list = as.character(rownames(degs))
  return(deg_list)
}

#3. convert essemble to entrez 
essemble_to_entrez = function (deg_list){
  entrez_map = as.data.frame(unlist(as.list(org.Mm.egENSEMBL2EG)))
  ind<-match(deg_list, rownames(entrez_map))
  deg_eg = entrez_map[ind,]
  return(deg_eg)
}

#get result
data = preprocess()
vsd <- getVarianceStabilizedData(data)
pdf("correlation.pdf", width=12, height=8)
pheatmap(cor(vsd),cexCol=0.75,cexRow=0.75)
dev.off()





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





data2 = preprocess_2()
contrast = c("Time", "inflam", "control")
res_dds_main=na.omit(as.data.frame(results(data2, contrast = contrast, cooksCutoff=FALSE)))
deg_inflam_VS_ctrl = res_dds_main[,c("log2FoldChange", "padj")]
contrast = c("Time", "tumor", "control")
res_dds_main=na.omit(as.data.frame(results(data2, contrast = contrast, cooksCutoff=FALSE)))
deg_tumor_VS_ctrl = res_dds_main[,c("log2FoldChange", "padj")]
contrast = c("Time", "tumor", "inflam")
res_dds_main=na.omit(as.data.frame(results(data2, contrast = contrast, cooksCutoff=FALSE)))
deg_tumor_VS_inflam = res_dds_main[,c("log2FoldChange", "padj")]



#--------------------------construct data------------------------------------------------------------------------------------------------
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


write.table(merge_data,"all_diff_data.txt", sep="\t", quote=F, row.names=T, col.names = T, na="NA", eol="\n")








degs = res_dds_main[res_dds_main$log2FoldChange > 1 & res_dds_main$padj < 0.05,]


#---------------------------------------------画图---------------------------------------------------------------
data = preprocess()


degfunc = function (group="", data, contrast=""){
  if (group != ""){
    res_dds_main=na.omit(as.data.frame(results(data, name=group, cooksCutoff=FALSE)))
    print("111111")
  }
  if (contrast != ""){
    res_dds_main=na.omit(as.data.frame(results(data, contrast = contrast, cooksCutoff=FALSE)))
    print ("222222")
  }
  degs = res_dds_main[res_dds_main$log2FoldChange < -1 & res_dds_main$padj < 0.05,]
  deg_list = as.character(rownames(degs))
  return(deg_list)
}


deg_list = degfunc(group = "Time_2Week_vs_control", data=data) 
length(deg_list)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("2wkVSctrl.down.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

deg_list = degfunc(contrast = c("Time", "4Week", "control"), data=data)
length(deg_list)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("4wkVSctrl.down.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

deg_list = degfunc(contrast = c("Time", "7Week", "control"), data=data)
length(deg_list)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("7wkVSctrl.down.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


deg_list = degfunc(contrast = c("Time", "10Week", "control"), data=data)
length(deg_list)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("10wkVSctrl.down.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()





data2 = preprocess_2()

deg_list = degfunc(contrast = c("Time", "inflam", "control"), data=data2)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("inflamVSctrl.down.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

deg_list = degfunc(contrast = c("Time", "tumor", "control"), data=data2)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("tumorVSctrl.down.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

deg_list = degfunc(contrast = c("Time", "tumor", "inflam"), data=data2)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("tumorVSinflam.down.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()



degfunc = function (group="", data, contrast=""){
  if (group != ""){
    res_dds_main=na.omit(as.data.frame(results(data, name=group, cooksCutoff=FALSE)))
    print("111111")
  }
  if (contrast != ""){
    res_dds_main=na.omit(as.data.frame(results(data, contrast = contrast, cooksCutoff=FALSE)))
    print ("222222")
  }
  degs = res_dds_main[res_dds_main$log2FoldChange > 1 & res_dds_main$padj < 0.05,]
  deg_list = as.character(rownames(degs))
  return(deg_list)
}


deg_list = degfunc(contrast = c("Time", "inflam", "control"), data=data2)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("inflamVSctrl.up.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

deg_list = degfunc(contrast = c("Time", "tumor", "control"), data=data2)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("tumorVSctrl.up.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

deg_list = degfunc(contrast = c("Time", "tumor", "inflam"), data=data2)
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("tumorVSinflam.up.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


#------------------------------------------------finish----------------------------------------------------

deg_list = degfunc(group = "Time_4Week_vs_control", data=data)









entrez_genes = essemble_to_entrez(deg_list = deg_list)
gene_list = as.vector(na.omit(entrez_genes))
write.table(gene_list, "/home/zhihl/Project/CRC/gene_list4.txt", quote =F , sep="\t", row.names=F)
#result in more than two group
#results(data, contrast = c("Time", "4Week", "control"))

#4. plot vocano 
vacano_plot = function(data){
  data1 = as.data.frame(results(data,name="Time_2Week_vs_control", cooksCutoff=FALSE))
  data1$change <-  as.factor(ifelse(data1$padj < 0.05 & abs(data1$log2FoldChange) > 1,ifelse(data1$log2FoldChange > 1,'UP','DOWN'),'NOT'))
  ggplot(data = data1, aes(x = -log10(padj), y = log2FoldChange, color = change)) +
    geom_point(alpha=0.8, size = 1) +
    theme_bw(base_size = 15) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    ) +
    scale_color_manual(name = "", values = c("deeppink1", "lightblue2", "snow2"), limits = c("UP", "DOWN", "NOT"))

}

vacano_plot(data)


#correlation plot
all_FPKM <- read.table("/home/zhihl/Project/CRC/genes.fpkm_table",header = T,sep = "\t")
track_id = all_FPKM[1]
control = all_FPKM[1:4]
control$mean_control = rowMeans(control[2:4])
control_mean_fpkm<- subset(control,select=c("tracking_id","mean_control"))


week2 = cbind(track_id, all_FPKM[5:7])
week2$mean_week2 = rowMeans(week2[2:4])
week2_mean_fpkm = subset(week2, select=c("tracking_id","mean_week2"))

week4 = cbind(track_id, all_FPKM[8:10])
week4$mean_week4 = rowMeans(week4[2:4])
week4_mean_fpkm = subset(week4, select=c("tracking_id","mean_week4"))

week7 = cbind(track_id, all_FPKM[11:13])
week7$mean_week7 = rowMeans(week7[2:4])
week7_mean_fpkm = subset(week7, select=c("tracking_id","mean_week7"))

week10 = cbind(track_id, all_FPKM[14:16])
week10$mean_week10 = rowMeans(week10[2:4])
week10_mean_fpkm = subset(week10, select=c("tracking_id","mean_week10"))

#control VS week2
expression_diff = function (data1, data2){
control_mean_fpkm = data1
week2_mean_fpkm = data2
control_week2<- merge(control_mean_fpkm,week2_mean_fpkm,by="tracking_id")
control_week2$log2_mean_control = log2(control_week2$mean_control + 1)
control_week2$log2_mean_week2 = log2(control_week2$mean_week7 + 1)
control_week2$FC1<- control_week2$mean_control/control_week2$mean_week7
control_week2$FC2<- control_week2$mean_week7/control_week2$mean_control

figure1 = control_week2
figure1$change = "NO"
control_hight = which(figure1$mean_control>5 & figure1$FC1 > 3)
figure1[control_hight,8]<- "control_hight"
week2_hight = which(figure1$mean_week2>5 & figure1$FC2 > 3)
figure1[week2_hight,8]<- "week7_hight"

#plot
ggplot(figure1,aes(x=log2_mean_control,y=log2_mean_week2,color = change))+ geom_point()+scale_color_manual(values = c("red", "lightgrey", "navy"))
control_high<- figure1[figure1$mean_control>5&figure1$FC1>3,]
week2_high<- figure1[figure1$mean_week7>5&figure1$FC2>3,]
data_different<- subset(figure1,figure1$mean_control>5|figure1$mean_week2>5)
#control_high=246,week2_hight=521,data_different=10080,
#线性
lm.data<- lm(figure1$log2_mean_week2~figure1$log2_mean_control)
summary(lm.data)
#R2 = 0.9473, beta:1.038596
}

expression_diff(data1=control_mean_fpkm, data2=week2_mean_fpkm)
expression_diff(data1=control_mean_fpkm, data2=week7_mean_fpkm)


#rpm
rpm = function (){
  setwd("/home/zhihl/Project/CRC")
  total_length =c(148013706, 79941339, 165711979, 134292652, 165600600, 93766819, 169674625, 97577473, 152838333, 110293532, 136734657, 97028485, 107600837, 98208097, 109510910)
  d.raw <- read.delim("rpm_count_table.txt",sep = "\t")
  d <- d.raw[rowSums(d.raw>3) > 2,]
  readcounts = subset(d, select=c("X823_control_1.txt", "X824_control_2.txt", "X825_control_3.txt", "X826_2week_1.txt", "X828_2week_2.txt", "X829_2week_3.txt"))
  names(readcounts)=c("control_1","control_2","control_3","week2_1","week2_2","week2_3")
  n = 0
  for (i in colnames(readcounts)){
    n = n + 1
    readcounts[,paste("rpm",i,sep="_")]=readcounts[,i]/total_length[n]*1000000
  }
  readcounts$rmp_control_mean<- 1/3*(readcounts[,7] + readcounts[,8] + readcounts[,9])
  readcounts$rmp_week2_mean<- 1/3*(readcounts[,10] + readcounts[,11] + readcounts[,12])
  #readcounts$p= apply(readcounts[,c("rmp_week2_mean","rmp_control_mean")],1,function(a){fisher.test(matrix(c(a,1000000 ,1000000),nrow = 2))$p})
  #readcounts$p= apply(readcounts[,c("rmp_week2_mean","rmp_control_mean")],1,function(a){fisher.test(matrix(c(a,1000000 -a[1] ,1000000 - a[2]),nrow = 2))$p})
  readcounts$FC<- (readcounts[,14]/readcounts[,13])  
  RRR_r=na.omit(readcounts[readcounts$FC<=1/5 ,])
  PRR_r=na.omit(readcounts[readcounts$FC>1/5 & readcounts$FC < 0.5 ,])
  FRR=na.omit(readcounts[readcounts$FC<=2 & readcounts$FC >=0.5,])
  PRR=na.omit(readcounts[readcounts$FC>2 & readcounts$FC<=5,])
  RRR=na.omit(readcounts[readcounts$FC>5,])

  region_supress = RRR_r[RRR_r$rmp_week2_mean > 10 , ]
  region_1 = t(as.data.frame(apply(as.data.frame(rownames(region_supress)),1, function(a){strsplit(a, "_")})))
  write.table(region_1, file="region_supress.txt", quote = F, sep="\t", row.names = FALSE, col.names = FALSE)
  
  region_less_supress = PRR_r[PRR_r$rmp_week2_mean > 10 , ]
  region_2 = t(as.data.frame(apply(as.data.frame(rownames(region_less_supress)),1, function(a){strsplit(a, "_")})))
  write.table(region_2, file="region_less_supress.txt", quote = F, sep="\t", row.names = FALSE, col.names = FALSE)
  
  region_active = RRR[RRR$rmp_week2_mean > 10 , ]
  region_3 = t(as.data.frame(apply(as.data.frame(rownames(region_active)),1, function(a){strsplit(a, "_")})))
  write.table(region_3, file="region_active.txt", quote = F, sep="\t", row.names = FALSE, col.names = FALSE)
  
  region_less_active = PRR[PRR$rmp_week2_mean > 10 , ]
  region_4 = t(as.data.frame(apply(as.data.frame(rownames(region_less_active)),1, function(a){strsplit(a, "_")})))
  write.table(region_4, file="region_less_active.txt", quote = F, sep="\t", row.names = FALSE, col.names = FALSE)
  
  region_no_change = FRR[FRR$rmp_week2_mean > 10 , ]
  region_5 = t(as.data.frame(apply(as.data.frame(rownames(region_no_change)),1, function(a){strsplit(a, "_")})))
  write.table(region_5, file="region_no_change.txt", quote = F, sep="\t", row.names = FALSE, col.names = FALSE)
  
}


compare_rpm = function(){
  setwd("/home/zhihl/Project/CRC/active")
  total_length =c(148013706, 79941339, 165711979, 134292652, 165600600, 93766819, 169674625, 97577473, 152838333, 110293532, 136734657, 97028485, 107600837, 98208097, 109510910)
  d.raw <- read.delim("rpm_count_table.txt",sep = "\t")
  readcounts = subset(d.raw, select=c("X826_2week_1.txt", "X828_2week_2.txt", "X829_2week_3.txt", "X830_4week_1.txt", "X831_4week_2.txt", "X832_4week_3.txt" ))
  names(readcounts)=c("week2_1","week2_2","week2_3", "week4_1","week4_2","week4_3")
  n = 0
  for (i in colnames(readcounts)){
    n = n + 1
    readcounts[,paste("rpm",i,sep="_")]=readcounts[,i]/total_length[n + 3]*1000000
  }
  readcounts$rmp_week2_mean<- 1/3*(readcounts[,7] + readcounts[,8] + readcounts[,9])
  readcounts$rmp_week4_mean<- 1/3*(readcounts[,10] + readcounts[,11] + readcounts[,12])
  readcounts$FC<- (readcounts[,14]/readcounts[,13])
}

#====================================================================
res <- results(d.deseq,name="Time_10Week_vs_control")
res2<-res[!is.na(res$padj) & res$padj<=0.01,]
deg<-as.character(rownames(res2))
library("org.Mm.eg.db")
ls("package:org.Mm.eg.db")
xx <- as.list(org.Mm.egENSEMBL2EG)
xxd <- as.data.frame(unlist(xx))
ind<-match(deg, rownames(xxd))

degeg<-xxd[ind,]
degeg

xx <- as.list(org.Mm.egGO)
ind<-match(degeg, names(xx))
degeggo<-xx[ind]

degeggo<-select(org.Mm.eg.db, as.character(degeg),
                "GO", keytype = "ENTREZID")


library("clusterProfiler")
ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("2wkVSctrl.up.2f.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

sig <- res[which(res$padj < 0.01),]
sig.deseq <- rownames(sig)
length(sig.deseq)

res <- results(d.deseq,name="Time_7Week_vs_control")
sig <- res[which(res$padj < 0.01),]
sig.deseq <- rownames(sig)
length(sig.deseq)

res <- results(d.deseq,name="Time_4Week_vs_control")
sig <- res[which(res$padj < 0.01),]
sig.deseq <- rownames(sig)
length(sig.deseq)

res <- results(d.deseq,name="Time_2Week_vs_control")
sig <- res[which(res$padj < 0.01),]
sig.deseq <- rownames(sig)
length(sig.deseq)



vsd <- getVarianceStabilizedData(d.deseq)
heatmap(cor(vsd),cexCol=0.75,cexRow=0.75)
pheatmap(cor(vsd),cexCol=0.75,cexRow=0.75)


pr <- prcomp(t(vsd))
plot(pr$x, main="PC plot",cex=0.1)
text(pr$x[,1],pr$x[,2],labels=colnames(vsd),cex=0.7)
