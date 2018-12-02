#set work dir
setwd("/home/zhihl/Project/CRC")
#include package, thses packages must be installed before
library(DESeq2)
library("org.Mm.eg.db")
library(ggplot2)


#1. process the read count data, and run deseq2
preprocess = function (){
  d.raw <- read.delim("count_table.txt",sep = "\t")
  d <- d.raw[rowSums(d.raw>3) > 2,]
  grp <- c(rep("control",3),rep("2Week",3),rep("4Week",3),rep("7Week",3),rep("10Week",3))
  cData <- data.frame(Time = factor(grp, levels = c("control", "2Week", "4Week", "7Week", "10Week")))
  rownames(cData) <- colnames(d)
  d.deseq <- DESeqDataSetFromMatrix(countData = d, colData = cData,design = ~Time)
  d.deseq <- DESeq(d.deseq)
  return(d.deseq)
}

#normalize
b=counts(data,normalized=T)

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
  degs = res_dds_main[res_dds_main$log2FoldChange > 2 & res_dds_main$padj < 0.05,]
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

data = preprocess()
#deg_list = degfunc(group = "Time_2Week_vs_control", data=data) 
#deg_list = degfunc(contrast = c("Time", "4Week", "control"), data=data)
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


control_week2<- merge(control_mean_fpkm,week2_mean_fpkm,by="tracking_id")
control_week2$log2_mean_control = log2(control_week2$mean_control + 1)
control_week2$log2_mean_week2 = log2(control_week2$mean_week2 + 1)
control_week2$FC1<- control_week2$mean_control/control_week2$mean_week2
control_week2$FC2<- control_week2$mean_week2/control_week2$mean_control

figure1 = control_week2
figure1$change = "NO"
control_hight = which(figure1$mean_control>5 & figure1$FC1 > 3)
figure1[control_hight,8]<- "control_hight"
week2_hight = which(figure1$mean_week2>5 & figure1$FC2 > 3)
figure1[week2_hight,8]<- "week2_hight"

#plot
ggplot(figure1,aes(x=log2_mean_control,y=log2_mean_week2,color = change))+ geom_point()+scale_color_manual(values = c("red", "lightgrey", "navy"))
control_high<- figure1[figure1$mean_control>5&figure1$FC1>3,]
week2_high<- figure1[figure1$mean_week2>5&figure1$FC2>3,]
data_different<- subset(figure1,figure1$mean_control>5|figure1$mean_week2>5)
#control_high=246,week2_hight=521,data_different=10080,
#线性
lm.data<- lm(figure1$log2_mean_week2~figure1$log2_mean_control)
summary(lm.data)
#R2 = 0.9473, beta:1.038596




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
ego <- enrichGO(gene          = deg,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

dotplot(ego, showCategory=15)

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
