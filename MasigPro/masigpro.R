library("maSigPro")
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
#our data------------------------------------------------------------
Time = c(rep(0, 18), rep(2, 18), rep(4, 18), rep(7, 18), rep(10, 18)) 
Replicates = rep(1:30, each=3)
Control = c(rep(1,18), rep(0, 72))
K27ac = rep(c(rep(1,3), rep(0,15)), 5)
K27me3 = rep(c(rep(0,3)), c(rep(1,3), rep(0,12)), 5)
K4me1 = rep(c(rep(0,6)), c(rep(1,3), rep(0,9)), 5)
K4me3 = rep(c(rep(0,9)), c(rep(1,3), rep(0,6)), 5)
K9me2 = rep(c(rep(0,12)), c(rep(1,3), rep(0,3)), 5)
K9me3 = rep(c(rep(0,15)), c(rep(1,3)), 5)

CRC.design = cbind(Time,Replicates, Control, K27ac,K27me3, K4me1, K4me3, K9me2, K9me3)
rownames(CRC.design) <- paste("Array", c(1:90), sep = "")
d <- make.design.matrix(CRC.design, degree = 5)
d
