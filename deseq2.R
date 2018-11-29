setwd("/home/zhihl/Project/CRC")

d.raw <- read.delim("count_table.txt",sep = "\t")
d <- d.raw[rowSums(d.raw>3) > 2,]
library(DESeq2)
grp <- c(rep("control",3),rep("2Week",3),rep("4Week",3),rep("7Week",3),rep("10Week",3))
cData <- data.frame(Time = factor(grp, levels = c("control", "2Week", "4Week", "7Week", "10Week")))
rownames(cData) <- colnames(d)
d.deseq <- DESeqDataSetFromMatrix(countData = d, colData = cData,design = ~Time)
d.deseq <- DESeq(d.deseq)

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


##########################################################

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

