library("maSigPro")
if (FALSE){
data(edesign.abiotic)
data(data.abiotic)
design <- make.design.matrix(edesign.abiotic, degree = 2)
fit <- p.vector(data.abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)
fit$i
fit$SELEC
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")
names(sigs)
names(sigs$sig.genes)
names(sigs$sig.genes$ColdvsControl)
suma2Venn(sigs$summary[, c(2:4)])
suma2Venn(sigs$summary[, c(1:4)])
sigs$sig.genes$SaltvsControl$g
pdf("hello.pdf")
see.genes(sigs$sig.genes$ColdvsControl, show.fit = F, dis =design$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 9, newX11 = FALSE)

dev.off()
STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]
PlotGroups (STMDE66, edesign = edesign.abiotic)
PlotGroups (STMDE66, edesign = edesign.abiotic, show.fit = T,  dis = design$dis, groups.vector = design$groups.vector)


data(edesignCT)

data(NBdata)
data(NBdesign)
d <- make.design.matrix(NBdesign)
library(MASS)
NBp <- p.vector(NBdata, d, counts=TRUE)
NBt <- T.fit(NBp)
get<-get.siggenes(NBt, vars="all")
get$summary
pdf("hello_rna.pdf")
see.genes(get$sig.genes, k = 4, newX11 = FALSE)
dev.off()
NBp <- p.vector(NBdata, d, counts=TRUE, theta=5)
NBp <- p.vector(NBdata, d, family=poisson() )


data(edesignCT)


ss.GENE <- function(n, r, var11 = 0.01, var12 = 0.02, var13 = 0.02,
                     a1 = 0, a2 = 0, a3 = 0) {
  tc.dat <- NULL
  for (i in 1:n) {
    gene <- c(rnorm(r, a1, var11), rnorm(r, a1, var12),
              rnorm(r, a3, var13))
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat }

flat <-ss.GENE(n = 850, r = 30) # flat
induc <- ss.GENE(n = 50, r = 30, a1 = 0, a2 = 0.2, a3 = 0.6) # induction
sat <- ss.GENE(n = 50, r = 30, a1 = 0, a2 = 1, a3 = 1.1) # saturation
ord <- ss.GENE(n = 50, r = 30, a1 = -0.8, a2 = -1, a3 = -1.3) # intercept
ss.DATA <- rbind(flat, induc,sat,ord)
rownames(ss.DATA) <- paste("feature", c(1:1000), sep = "")
colnames(ss.DATA) <- paste("Array", c(1:90), sep = "")
}



#our data------------------------------------------------------------
#construct matrix
Time = c(rep(0, 18), rep(2, 18), rep(4, 18), rep(7, 18), rep(10, 18)) 
Replicates = rep(1:30, each=3)
Control = c(rep(1,18), rep(0, 72))
K27ac = rep(c(rep(1,3), rep(0,15)), 5)
K27me3 = rep(c(rep(0,3), rep(1,3), rep(0,12)), 5)
K4me1 = rep(c(rep(0,6), rep(1,3), rep(0,9)), 5)
K4me3 = rep(c(rep(0,9), rep(1,3), rep(0,6)), 5)
K9me2 = rep(c(rep(0,12), rep(1,3), rep(0,3)), 5)
K9me3 = rep(c(rep(0,15), rep(1,3)), 5)

CRC.design = cbind(Time,Replicates, Control, K27ac,K27me3, K4me1, K4me3, K9me2, K9me3)
#rownames(CRC.design) <- paste("Array", c(1:90), sep = "")
sample_vector = c()
for (we in c("ctrl", "2weeks", "4weeks", "7weeks", "10weeks")){
  for (mark in c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me3", "H3K9me2", "H3K9me3")){
    for ( rep in c("1", "2" ,"3")){
      sample = paste(paste(we, rep, sep="_"), mark, sep="_")
      sample_vector = c(sample_vector, sample)
    }
  }
}
rownames(CRC.design) = sample_vector
d <- make.design.matrix(CRC.design, degree = 3)
d



#training data
df = read.delim("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/expression_10000_ensembl.txt", sep=",", header=T, row.names=1, check.names=FALSE)
df = df[1:700, ]
df = data.matrix(df)
df = scale(df, center = TRUE, scale = TRUE)
#library(MASS)
#NBp <- p.vector(df, d, counts=TRUE)
#NBt <- T.fit(NBp)
fit <- p.vector(df, d, Q = 0.05, MT.adjust = "BH", min.obs = 20)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)

get<-get.siggenes(tstep, vars="all")
pdf("hello_rna_version2.pdf")
cluster_result = see.genes(get$sig.genes, k = 9, newX11 = FALSE)
dev.off()

#save.image(file = "/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/test.RData")



#----
if (FALSE){
fit <- p.vector(df, d, Q = 0.05, MT.adjust = "BH", min.obs = 20)
tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
#save(tstep, file = "/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/masigpro_all_degree3.RData")
get<-get.siggenes(tstep, vars="all")
pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/cluster_result.pdf")
see.genes(get$sig.genes, k = 20, newX11 = FALSE)
dev.off()

load("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/masigpro_all_degree3.RData")
get<-get.siggenes(tstep, vars="groups")
suma2Venn(get$summary[, c(2:7)])
pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/H3k27ac_result.pdf")
cluster_result = see.genes(get$sig.genes$K27acvsControl, dis =d$dis,
          cluster.method="hclust" ,cluster.data = 1, k = 9, newX11 = FALSE)
dev.off()


get$sig.genes$K27me3vsControl$g
get$sig.genes$K27acvsControl$g
get$sig.genes$K4me3vsControl$g
get$sig.genes$K4me1vsControl$g
get$sig.genes$K9me2vsControl$g
get$sig.genes$K9me3vsControl$g
H3k27ac_genes_list = rownames(get$sig.genes$K27acvsControl$sig.profiles)
H3k27me3_genes_list = rownames(get$sig.genes$K27me3vsControl$sig.profiles)
H3k4me3_genes_list = rownames(get$sig.genes$K4me3vsControl$sig.profiles)
H3k9me3_genes_list = rownames(get$sig.genes$K9me3vsControl$sig.profiles)

length(H3k27me3_genes_list)
length(H3k4me3_genes_list)
length(H3k9me3_genes_list)

length(intersect(H3k27me3_genes_list, total_H3K27me3))
length(intersect(H3k9me3_genes_list, total_H3K9me3))
length(intersect(H3k4me3_genes_list, total_H3K4me3))
}
#-------------------------------------
if (FALSE){
load("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/masigpro_all_degree3.RData")
get<-get.siggenes(tstep, vars="all")

pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/cluster_result.pdf")
cluster_result = see.genes(get$sig.genes, k = 20, newX11 = FALSE)
dev.off()

clu8 = cluster_result$cut[cluster_result$cut == 8]
clu9 = cluster_result$cut[cluster_result$cut == 9]
clu15 = cluster_result$cut[cluster_result$cut == 15]
deg_list = c(names(clu8), names(clu9), names(clu15))

ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/cluster8_9_15.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


clu1 = cluster_result$cut[cluster_result$cut == 1]
clu5 = cluster_result$cut[cluster_result$cut == 5]
clu13 = cluster_result$cut[cluster_result$cut == 13]
clu18 = cluster_result$cut[cluster_result$cut == 18]
deg_list = c(names(clu1), names(clu5), names(clu13), names(clu18))

ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/cluster1_5_13_18.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


clu3 = cluster_result$cut[cluster_result$cut == 3]
clu4 = cluster_result$cut[cluster_result$cut == 4]
clu6 = cluster_result$cut[cluster_result$cut == 6]
clu17 = cluster_result$cut[cluster_result$cut == 17]
clu19 = cluster_result$cut[cluster_result$cut == 19]
deg_list = c(names(clu3), names(clu4), names(clu6), names(clu17), names(19))

ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/cluster3_4_6_17_19.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()



clu7 = cluster_result$cut[cluster_result$cut == 7]
clu10 = cluster_result$cut[cluster_result$cut == 10]
clu11 = cluster_result$cut[cluster_result$cut == 11]
clu20 = cluster_result$cut[cluster_result$cut == 20]
deg_list = c(names(clu7), names(clu10), names(clu11), names(clu20))

ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/cluster7_10_11_20.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


clu2 = cluster_result$cut[cluster_result$cut == 2]
clu12 = cluster_result$cut[cluster_result$cut == 12]
clu14 = cluster_result$cut[cluster_result$cut == 14]
clu16 = cluster_result$cut[cluster_result$cut == 16]
deg_list = c(names(clu2), names(clu12), names(clu14), names(clu16))

ego <- enrichGO(gene          = deg_list,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/cluster2_12_14_16.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()




ego <- enrichGO(gene          = total_H3K27ac_cancer,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/total_H3K27ac_cancer.pdf", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()


ego <- enrichGO(gene          = total_H3K27ac_colits,
                keyType       = "ENSEMBL",
                OrgDb         = "org.Mm.eg.db",
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

pdf("/home/zhihl/Project/CRC/Chip_analysis/MaSigPro/total_H3K27ac_colits", width=10, height=10)
dotplot(ego, showCategory=15)
dev.off()

}