#### 450k MethyLumi - Lumi - BMIQ - ComBat ####
# Endurance exercise: Before, After
setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/GSE60655_End')
# dataDir <- ('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/DataDir/')

##### Libraries #####
suppressMessages(
    suppressWarnings({
        library(GEOquery)
        
        library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        library(IlluminaHumanMethylation450kmanifest)
        library(lumi)  # MethyLumi
        library(wateRmelon)  # BMIQ
        library(sva)  # ComBat
        library(limma)
        library(missMethyl)
        library(RColorBrewer)
        library(stringr)
        library(dplyr)
        library(forcats)
        library(EnhancedVolcano)
    })
)

pal <- brewer.pal(8,"Dark2")  # for plotting


##### Obtaining the Data #####
# library(GEOquery)
# getGEOSuppFiles("GSE60655", makeDirectory = F) # 450k
# untar("GSE60655_RAW.tar", exdir = "idat")
# head(list.files("idat", pattern = "idat"))

# idatFiles <- list.files("idat", pattern = "idat.gz$", full = TRUE)
# sapply(idatFiles, gunzip, overwrite = TRUE)

# files <- list.files("idat", pattern = ".idat", full = TRUE); length(files)
# str_view(files, "GSM[0-9]+_?")
# files.renamed <- gsub("GSM[0-9]+_?", "", files)
# file.rename(files, files.renamed)

idatFiles <- list.files("idat", pattern = ".idat", full.names = TRUE)
barcodes <- str_extract(idatFiles, "\\d+_R\\d+C\\d+"); length(barcode)  # 72
mlset <- methylumIDAT(barcodes = barcode, idatPath = "idat"); dim(mlset)  # MethyLumiSet, 485577  36
mlset.raw <- mlset; mlset # mlset <- mlset.raw
# Object Information:
# MethyLumiSet (storageMode: lockedEnvironment)
# assayData: 485577 features, 36 samples 
# element names: betas, methylated, methylated.OOB, pvals, unmethylated, unmethylated.OOB 
# protocolData: none
# phenoData
# sampleNames: 5684819044_R01C01 5684819044_R01C02 ... 6164621126_R06C02 (36 total)
# varLabels: barcode
# varMetadata: labelDescription
# featureData
# featureNames: cg00000029 cg00000108 ... rs9839873 (485577 total)
# fvarLabels: Probe_ID DESIGN COLOR_CHANNEL
# fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: IlluminaHumanMethylation450k 
# Major Operation History:
#     submitted            finished                                             command
# 1 2021-08-31 15:49:46 2021-08-31 15:51:35 methylumIDAT(barcodes = barcode, idatPath = "idat")
# 2 2021-08-31 15:51:40 2021-08-31 15:51:46                          Subset of 485577 features.

write.table(featureNames(mlset), "GSE60655.org.probes.txt", quote = F,
            row.names = F, col.names = F)
# test <- read.table("GSE60655.org.probes.txt", header = F)
# head(test,5)
# nrow(test)


mlset@assayData$betas[1:3,1:5]
#            5684819044_R01C01 5684819044_R01C02 5684819044_R02C01 5684819044_R02C02 5684819044_R03C01
# cg00000029         0.4323047         0.3802386         0.4390304         0.4111430         0.3364002
# cg00000108         0.9205592         0.8975848         0.8969463         0.9230938         0.9157878
# cg00000109         0.8125258         0.7672213         0.8282256         0.7944335         0.7895413

mlset@assayData$pvals[1:3,1:5]
#            5684819044_R01C01 5684819044_R01C02 5684819044_R02C01 5684819044_R02C02 5684819044_R03C01
# cg00000029                 0                 0                 0                 0                 0
# cg00000108                 0                 0                 0                 0                 0
# cg00000109                 0                 0                 0                 0                 0


source("pdata.sort.R"); dim(pD)  # 14  11
head(pD, 3)



#### Adjusting the order of rownames(pD) & featureNames(mlset) ####
table(seq(1,nrow(pD)) == match(pD$barcode, sampleNames(mlset)))
# FALSE  TRUE 
# 11     3


sampleNames(mlset)
# [1] "5684819044_R01C01" "5684819044_R01C02" "5684819044_R02C01" "5684819044_R02C02" "5684819044_R03C01" "5684819044_R03C02" "5684819044_R04C01"
# [8] "5684819044_R04C02" "5684819044_R05C01" "5684819044_R05C02" "5684819044_R06C01" "5684819044_R06C02" "6164621122_R01C01" "6164621122_R01C02"
# [15] "6164621122_R02C01" "6164621122_R02C02" "6164621122_R03C01" "6164621122_R03C02" "6164621122_R04C01" "6164621122_R04C02" "6164621122_R05C01"
# [22] "6164621122_R05C02" "6164621122_R06C01" "6164621122_R06C02" "6164621126_R01C01" "6164621126_R01C02" "6164621126_R02C01" "6164621126_R02C02"
# [29] "6164621126_R03C01" "6164621126_R03C02" "6164621126_R04C01" "6164621126_R04C02" "6164621126_R05C01" "6164621126_R05C02" "6164621126_R06C01"
# [36] "6164621126_R06C02"


pD$barcode
# [1] "5684819044_R01C01" "5684819044_R02C01" "5684819044_R03C01" "5684819044_R04C01" "5684819044_R05C01" "5684819044_R06C01" "5684819044_R01C02"
# [8] "5684819044_R02C02" "5684819044_R03C02" "5684819044_R04C02" "5684819044_R05C02" "5684819044_R06C02" "6164621122_R01C01" "6164621122_R02C01"


sampleNames(mlset)[match(pD$barcode, sampleNames(mlset))]
# [1] "5684819044_R01C01" "5684819044_R02C01" "5684819044_R03C01" "5684819044_R04C01" "5684819044_R05C01" "5684819044_R06C01" "5684819044_R01C02"
# [8] "5684819044_R02C02" "5684819044_R03C02" "5684819044_R04C02" "5684819044_R05C02" "5684819044_R06C02" "6164621122_R01C01" "6164621122_R02C01"


mlset <- mlset[,match(pD$barcode, sampleNames(mlset))]
sampleNames(mlset); length(sampleNames(mlset))  # 14
# [1] "5684819044_R01C01" "5684819044_R02C01" "5684819044_R03C01" "5684819044_R04C01" "5684819044_R05C01" "5684819044_R06C01" "5684819044_R01C02"
# [8] "5684819044_R02C02" "5684819044_R03C02" "5684819044_R04C02" "5684819044_R05C02" "5684819044_R06C02" "6164621122_R01C01" "6164621122_R02C01"


#### merging pData & mlset ####
pData(mlset) <- pD


#. renaming rownames(pD) & sampleNames(mlset) ----
rownames(pD) <- pD$ID
sampleNames(mlset) <- rownames(pD)
mlset
# MethyLumiSet (storageMode: lockedEnvironment)
# assayData: 485577 features, 14 samples 
# element names: betas, methylated, methylated.OOB, pvals, unmethylated, unmethylated.OOB 
# protocolData: none
# phenoData
# sampleNames: S1.before S1.after ... S7.after (14 total)
# varLabels: rep gender ... barcode (11 total)
# varMetadata: labelDescription
# featureData
# featureNames: cg00000029 cg00000108 ... rs9839873 (485577 total)
# fvarLabels: Probe_ID DESIGN COLOR_CHANNEL
# fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: IlluminaHumanMethylation450k 



#### removing 65 SNPs ####
length(grep("rs[0-9]+", featureNames(mlset)))  # 65
snps <- featureNames(mlset)[grep("rs[0-9]+", featureNames(mlset))]
keep <- !(featureNames(mlset) %in% snps); table(keep)
# keep
# FALSE   TRUE 
# 65 485512
mlset <- mlset[keep,]; dim(mlset)  # 485512  14


#### Removing probes on sex chromosomes (skipped) ####
#. Get 450k Annotation Data ----
# ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# keep <- !(featureNames(mlset) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]); table(keep)
# keep
# FALSE   TRUE 
# 11648 473864
# mlset <- mlset[keep,]; dim(mlset)  # 473864  34



#### detection p-val #### 
detP <- mlset@assayData$pvals; dim(detP)  # 473864  14
length(detP[detP > 0.01])  # 2062
max(colMeans(detP))  # 0.0001693573
detP[1:3, 1:5]
#            S1.before S1.after S4.before S4.after S5.before
# cg00000029         0        0         0        0         0
# cg00000108         0        0         0        0         0
# cg00000109         0        0         0        0         0


#. plot detP ----
# examine mean detection p-values across all samples to identify any failed samples
par(mfrow=c(1,1))
# All・Timepoint
barplot(colMeans(detP), col=pal[factor(pD$timepoint)], las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(pD$timepoint)), fill=pal, bg="white", cex = 0.8, bty = "n")


#. removing detection p value > 0.01 in > 5% of the samples ----
table(rowSums(detP < 0.01) > 0.95*ncol(mlset))  # 'FALSE' to be removed
# FALSE   TRUE 
#   978 484534

keep <- rowSums(detP < 0.01) > 0.95*ncol(mlset)
mlset <- mlset[keep,]; dim(mlset)  # 484534  14
detP <- detP[keep,]; dim(detP)  # 484534  14



#### removing outlier (b & M vals) ####
mVals.pre <- estimateM(as(mlset, "MethyLumiM"), returnType = "matrix")
bVals.pre <- mlset@assayData$betas

par(mfrow=c(1,2))
det.outlier.m.pre <- detectOutlier(mVals.pre, ifPlot = TRUE)  # none
det.outlier.b.pre <- detectOutlier(bVals.pre, ifPlot = TRUE)  # none



#### Color-bias adjustment ####
mlset.clradj <- adjColorBias.quantile(mlset)
par(mfrow = c(1,2))
plotColorBias1D(mlset)
plotColorBias1D(mlset.clradj)



#### BMIQ (probe type bias adjustment) ####
# normalization on β-values
# library(wateRmelon)
mlset.clradj.bmiq <- BMIQ(mlset.clradj)
par(mfrow = c(1,2))
plotDensity(mlset@assayData$betas, legend = FALSE, main = "raw", col=pal[factor(pD$ID)])
legend("topleft", legend=levels(factor(pD$ID)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")
plotDensity(mlset.clradj.bmiq@assayData$betas, legend = FALSE, main = "BMIQed", col=pal[factor(pD$ID)])
legend("topleft", legend=levels(factor(pD$ID)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")



#### Data exploration ####
#. MDS plot ----
# Multi-dimensional scaling (MDS) plots: to look at largest sources of variation
mVals <- estimateM(as(mlset.clradj.bmiq, "MethyLumiM"), returnType = "matrix")
bVals <- mlset.clradj.bmiq@assayData$betas

par(mfrow=c(1,4))
plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$timepoint)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$sample)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$sample)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$batch)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$batch)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$slide)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$slide)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")


#. Examine Higher dimensions ----
par(mfrow=c(1,4))
plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], dim=c(1,3))
legend("topleft", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], dim=c(2,3))
legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], dim=c(3,4))
legend("bottomright", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")


#. density plot: M & beta values ----
par(mfrow=c(1,2))
plotDensity(mVals, legend = FALSE, main = "M values", col=pal[factor(pD$timepoint)])
legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")
plotDensity(bVals, legend = FALSE, main = "Beta values", col=pal[factor(pD$timepoint)])
legend("topleft", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")



#### ComBat: removing batch effect ####
batch <- as.factor(pD$slide)
# mod <- model.matrix(~ as.factor(pD$timepoint) + as.factor(pD$sample), data = pD)
# mod0 <- model.matrix(~ as.factor(pD$sample), data = pD)
modcombat <- model.matrix(~ 1, data = pD)
combat.edata <- ComBat(dat = mVals, batch = batch, mod = modcombat, 
                       par.prior = TRUE, prior.plots = TRUE)
# combat.edata.npar <- ComBat(dat = mVals, batch = batch, mod = modcombat, 
#                        par.prior = FALSE, prior.plots = TRUE)
# pValuesComBat <- f.pvalue(combat.edata, mod, mod0)
# qValuesComBat <- p.adjust(pValuesComBat, method="BH")


#. MDS: visualizing Combat()ed data ----
par(mfrow=c(1,3))
plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$timepoint)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pD$sample)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$sample)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pD$batch)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$slide)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")



#### Probe-wise Differential Methylation Analysis ####
#. Create a Design Matrix ----
# ref: https://genomicsclass.github.io/book/pages/expressing_design_formula.html
timepoint <- factor(pD$timepoint)  # factor of interest
individual <- factor(pD$sample)  # adjustment
design <- model.matrix(~ 0 + timepoint + individual, data = pD)
colnames(design) <- c(levels(timepoint),levels(individual)[-1])
design[1:5, 1:5]

#. Fit the Linear Model ----
# lmFit: Linear Model For Series Of Arrays
fit <- lmFit(combat.edata, design)

#. Contrast Matrix ----
# for specific comparisons
contMatrix <- makeContrasts(after-before, levels=design)  # ref: before
# contMatrix

#. Fit the Contrasts ----
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

#. Numbers of DM CpGs at FDR < 0.05 ----
summary(decideTests(fit2))
summary(decideTests(fit2, p.value = 0.2))
summary(decideTests(fit2, adjust.method = "none", p.value = 0.05)) # p-value
# after - before (FDR < 0.05)
# Down                0
# NotSig         484534
# Up                  0

# after - before (FDR < 2.0)
# Down                0
# NotSig         484534
# Up                  0

# after - before
# Down            20589
# NotSig         436782
# Up              27163




#### Get the Table of Results ####
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(combat.edata), ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
# head(ann450kSub)

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub, sort.by = "p")
DMPs[1:3,c("logFC", "P.Value", "adj.P.Val")]
#                 logFC      P.Value adj.P.Val
# cg08642731  0.6338006 6.213410e-06 0.4997675
# cg02586025  0.5055961 2.447523e-05 0.4997675
# cg12337011 -0.5370560 3.103759e-05 0.4997675

DMPs %>% select(logFC, P.Value, adj.P.Val) %>% 
    filter(P.Value < 0.05 & abs(logFC) > 1) %>% 
    nrow()  # 3

summary(DMPs$logFC)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -1.222968 -0.097350  0.001802  0.004197  0.102325  1.169221


write.table(DMPs, file="DMPs.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv", sep=",", row.names=FALSE)
dmps.005 <- subset(DMPs, P.Value < 0.05)
write.table(dmps.005, file="DMPs.p005.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv", sep=",", row.names=FALSE)



#### Volcano plot ####
# library(EnhancedVolcano)
EnhancedVolcano(DMPs,
                lab = DMPs$Name,
                # selectLab = 'Peptide9A',
                x = 'logFC', FCcutoff = 0.5, #xlim = c(-2, 2),
                y = 'P.Value', pCutoff = 0.05, #ylim = c(0, 6),
                xlim = c(floor(min(DMPs$logFC)), ceiling(max(abs(DMPs$logFC)))),
                ylim = c(0, ceiling(max(-log10(DMPs$P.Value)))),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~italic(P)),  # ~adjusted
                pointSize = 2.0, 
                labSize = 3.0,
                col = c('black', 'red', 'orange', 'blue'),
                legendLabels = c('NS','Log2 FC','padj','padj & Log2 FC')) +

labs(title = "GSE60655 Endurance p < 0.05", titleLabSize = 10,
     subtitle = "pre vs post", caption = "")




#### Plot Top 4 Most Sig DM.CpGs ####
# cf. DMPs is already ordered by P.Value
# par(mfrow=c(2,2))
# sapply(rownames(dmps.005)[1:4], function(cpg){
#     plotCpg(bVals, cpg=cpg, pheno=pD$timepoint, ylab = "Beta values")
# })



#### Gene Ontology mlseting (gometh) ####
#. significant CpG sites at FDR < 0.05 ----
# library(missMethyl)
cpgs.sig.names <- dmps.005$Name; length(cpgs.sig.names)   # 47752
cpgs.pos.names <- dmps.005[dmps.005$logFC > 0,]$Name; length(cpgs.pos.names)  # 27163
cpgs.neg.names <- dmps.005[dmps.005$logFC < 0,]$Name; length(cpgs.neg.names)  # 20589
all.names <- DMPs$Name; length(all.names)  # 484534
gst <- gometh(sig.cpg=cpgs.sig.names, all.cpg=all.names, plot.bias=F)
gst.pos <- gometh(sig.cpg=cpgs.pos.names, all.cpg=all.names, plot.bias=F)
gst.neg <- gometh(sig.cpg=cpgs.neg.names, all.cpg=all.names, plot.bias=F)

#. categories ----
t10GO <- topGSA(gst, num=Inf)
t10GO <- t10GO[t10GO$P.DE < 0.05,]; length(t10GO$TERM)  # 1146
t10GO.pos <- topGSA(gst.pos, num=Inf)
t10GO.pos <- t10GO.pos[t10GO.pos$P.DE < 0.05,]; length(t10GO.pos$TERM)  # 482
t10GO.neg <- topGSA(gst.neg, num=Inf)
t10GO.neg <- t10GO.neg[t10GO.neg$P.DE < 0.05,]; length(t10GO.neg$TERM)  # 1336

t10GO$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO$TERM)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO$TERM)  # 0
t10GO.pos$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)  # 0
t10GO.neg$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.neg$TERM)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO.neg$TERM)  # 0


write.csv(t10GO, "go.p005.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv")
write.csv(t10GO.pos, "go.p005.pos.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv")
write.csv(t10GO.neg, "go.p005.neg.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv")
t10GO$TERM[1:20]
# [1] "anatomical structure morphogenesis"                                      "positive regulation of RNA biosynthetic process"                        
# [3] "positive regulation of nucleic acid-templated transcription"             "cell morphogenesis"                                                     
# [5] "cell morphogenesis involved in differentiation"                          "cellular component organization"                                        
# [7] "neuron projection morphogenesis"                                         "nucleoplasm"                                                            
# [9] "positive regulation of RNA metabolic process"                            "cellular component morphogenesis"                                       
# [11] "cellular component organization or biogenesis"                           "cell junction"                                                          
# [13] "positive regulation of transcription, DNA-templated"                     "animal organ morphogenesis"                                             
# [15] "plasma membrane bounded cell projection morphogenesis"                   "cell projection morphogenesis"                                          
# [17] "anatomical structure formation involved in morphogenesis"                "neuron differentiation"                                                 
# [19] "cell projection"                                                         "positive regulation of nucleobase-containing compound metabolic process"

