#### 450k MethyLumi - Lumi - BMIQ - ComBat ####
# chr remain, male only
# Endurance exercise: Before, After
setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/test/0902ref.before/GSE114763_Res')
# dataDir <- ('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/DataDir/')

##### Libraries #####
suppressMessages(
    suppressWarnings({
        library(GEOquery)
        library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
        library(IlluminaHumanMethylationEPICmanifest)
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

mlset <- readEPIC("idat", parallel = T, n = F); dim(mlset)  # MethyLumiSet, 866297  40
mlset.raw <- mlset
mlset
# MethyLumiSet (storageMode: lockedEnvironment)
# assayData: 866297 features, 40 samples 
# element names: betas, methylated, pvals, unmethylated 
# protocolData: none
# phenoData
# sampleNames: GSM3149860_201496860025_R01C01 GSM3149861_201496860025_R02C01
# ... GSM3149899_201496860128_R08C01 (40 total)
# varLabels: barcode
# varMetadata: labelDescription
# featureData
# featureNames: cg00000029 cg00000103 ... rs9839873 (866297 total)
# fvarLabels: Probe_ID DESIGN COLOR_CHANNEL
# fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: IlluminaHumanMethylationEpic 


#### Keeping the same probes with GSE60655 ####
# the probes are the original
keep.probes <- read.table("../GSE60655_End/GSE60655.org.probes.txt"); nrow(keep.probes)  # 485577
head(keep.probes)
keep.probe.names <- keep.probes$V1
table(featureNames(mlset) %in% keep.probe.names)
# FALSE   TRUE 
# 413406 452891

keep <- featureNames(mlset) %in% keep.probe.names
mlset <- mlset[keep,]; dim(mlset)  # 452891  40

mlset@assayData$betas[1:3,1:3]
#            GSM3149860_201496860025_R01C01 GSM3149861_201496860025_R02C01 GSM3149862_201496860025_R03C01
# cg00000029                      0.1434010                      0.1824713                      0.2594518
# cg00000109                      0.5752941                      0.6659634                      0.6598120
# cg00000165                      0.1014829                      0.1299783                      0.1196355

mlset@assayData$pvals[1:3,1:3]
#            GSM3149860_201496860025_R01C01 GSM3149861_201496860025_R02C01 GSM3149862_201496860025_R03C01
# cg00000029                              0                              0                              0
# cg00000109                              0                              0                              0
# cg00000165                              0                              0                              0


source("pdata.sort.R"); dim(pD)  # 16  9
head(pD, 3)
#             rep timepoint                                    file geo_accession sample          ID                        file.id             barcode        slide
# GSM3149860 rep1  baseline GSM3149860_201496860025_R01C01_Grn.idat    GSM3149860     S1 S1.baseline GSM3149860_201496860025_R01C01 201496860025_R01C01 201496860025
# GSM3149862 rep1   loading GSM3149862_201496860025_R03C01_Grn.idat    GSM3149862     S1  S1.loading GSM3149862_201496860025_R03C01 201496860025_R03C01 201496860025
# GSM3149865 rep1  baseline GSM3149865_201496860025_R06C01_Grn.idat    GSM3149865     S2 S2.baseline GSM3149865_201496860025_R06C01 201496860025_R06C01 201496860025



#### Adjusting the order of rownames(pD) & featureNames(mlset) ####
table(seq(1,nrow(pD)) == match(pD$file.id, sampleNames(mlset)))
# FALSE  TRUE 
# 15     1


sampleNames(mlset)
# [1] "GSM3149860_201496860025_R01C01" "GSM3149861_201496860025_R02C01" "GSM3149862_201496860025_R03C01"
# [4] "GSM3149863_201496860025_R04C01" "GSM3149864_201496860025_R05C01" "GSM3149865_201496860025_R06C01"
# [7] "GSM3149866_201496860025_R07C01" "GSM3149867_201496860025_R08C01" "GSM3149868_201465940053_R01C01"
# [10] "GSM3149869_201465940053_R02C01" "GSM3149870_201465940053_R03C01" "GSM3149871_201465940053_R04C01"
# [13] "GSM3149872_201465940053_R05C01" "GSM3149873_201465940053_R06C01" "GSM3149874_201465940053_R07C01"
# [16] "GSM3149875_201465940053_R08C01" "GSM3149876_201496850072_R01C01" "GSM3149877_201496850072_R02C01"
# [19] "GSM3149878_201496850072_R03C01" "GSM3149879_201496850072_R04C01" "GSM3149880_201496850072_R05C01"
# [22] "GSM3149881_201496850072_R06C01" "GSM3149882_201496850072_R07C01" "GSM3149883_201496850072_R08C01"
# [25] "GSM3149884_201496860106_R01C01" "GSM3149885_201496860106_R02C01" "GSM3149886_201496860106_R03C01"
# [28] "GSM3149887_201496860106_R04C01" "GSM3149888_201496860106_R05C01" "GSM3149889_201496860106_R06C01"
# [31] "GSM3149890_201496860106_R07C01" "GSM3149891_201496860106_R08C01" "GSM3149892_201496860128_R01C01"
# [34] "GSM3149893_201496860128_R02C01" "GSM3149894_201496860128_R03C01" "GSM3149895_201496860128_R04C01"
# [37] "GSM3149896_201496860128_R05C01" "GSM3149897_201496860128_R06C01" "GSM3149898_201496860128_R07C01"
# [40] "GSM3149899_201496860128_R08C01"


pD$file.id
# [1] "GSM3149860_201496860025_R01C01" "GSM3149862_201496860025_R03C01" "GSM3149865_201496860025_R06C01"
# [4] "GSM3149867_201496860025_R08C01" "GSM3149870_201465940053_R03C01" "GSM3149872_201465940053_R05C01"
# [7] "GSM3149875_201465940053_R08C01" "GSM3149877_201496850072_R02C01" "GSM3149880_201496850072_R05C01"
# [10] "GSM3149882_201496850072_R07C01" "GSM3149885_201496860106_R02C01" "GSM3149887_201496860106_R04C01"
# [13] "GSM3149890_201496860106_R07C01" "GSM3149892_201496860128_R01C01" "GSM3149895_201496860128_R04C01"
# [16] "GSM3149897_201496860128_R06C01"


sampleNames(mlset)[match(pD$file.id, sampleNames(mlset))]
# [1] "GSM3149860_201496860025_R01C01" "GSM3149862_201496860025_R03C01" "GSM3149865_201496860025_R06C01"
# [4] "GSM3149867_201496860025_R08C01" "GSM3149870_201465940053_R03C01" "GSM3149872_201465940053_R05C01"
# [7] "GSM3149875_201465940053_R08C01" "GSM3149877_201496850072_R02C01" "GSM3149880_201496850072_R05C01"
# [10] "GSM3149882_201496850072_R07C01" "GSM3149885_201496860106_R02C01" "GSM3149887_201496860106_R04C01"
# [13] "GSM3149890_201496860106_R07C01" "GSM3149892_201496860128_R01C01" "GSM3149895_201496860128_R04C01"
# [16] "GSM3149897_201496860128_R06C01"


mlset <- mlset[,match(pD$file.id, sampleNames(mlset))]
sampleNames(mlset); length(sampleNames(mlset))  # 16
# [1] "GSM3149860_201496860025_R01C01" "GSM3149862_201496860025_R03C01" "GSM3149865_201496860025_R06C01"
# [4] "GSM3149867_201496860025_R08C01" "GSM3149870_201465940053_R03C01" "GSM3149872_201465940053_R05C01"
# [7] "GSM3149875_201465940053_R08C01" "GSM3149877_201496850072_R02C01" "GSM3149880_201496850072_R05C01"
# [10] "GSM3149882_201496850072_R07C01" "GSM3149885_201496860106_R02C01" "GSM3149887_201496860106_R04C01"
# [13] "GSM3149890_201496860106_R07C01" "GSM3149892_201496860128_R01C01" "GSM3149895_201496860128_R04C01"
# [16] "GSM3149897_201496860128_R06C01"


#### merging pData & mlset ####
pData(mlset) <- pD


#. renaming rownames(pD) & sampleNames(mlset) ----
rownames(pD) <- pD$ID
sampleNames(mlset) <- rownames(pD)
mlset
# MethyLumiSet (storageMode: lockedEnvironment)
# assayData: 452891 features, 16 samples 
# element names: betas, methylated, pvals, unmethylated 
# protocolData: none
# phenoData
# sampleNames: S1.baseline S1.loading ... S8.loading (16 total)
# varLabels: rep timepoint ... file.id (9 total)
# varMetadata: labelDescription
# featureData
# featureNames: cg00000029 cg00000109 ... rs9839873 (452891 total)
# fvarLabels: Probe_ID DESIGN COLOR_CHANNEL
# fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: IlluminaHumanMethylationEpic



#### removing 65 SNPs ####
length(grep("rs[0-9]+", featureNames(mlset)))  # 59
snps <- featureNames(mlset)[grep("rs[0-9]+", featureNames(mlset))]
keep <- !(featureNames(mlset) %in% snps); table(keep)
# keep
# FALSE   TRUE 
# 59 452832
mlset <- mlset[keep,]; dim(mlset)  # 452832  16



#### Removing probes on sex chromosomes (skipped) ####
#. Get 450k Annotation Data ----
# ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# keep <- !(featureNames(mlset) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")]); table(keep)
# keep
# FALSE   TRUE 
# 11648 473864
# mlset <- mlset[keep,]; dim(mlset)  # 473864  34



#### detection p-val #### 
detP <- mlset@assayData$pvals; dim(detP)  # 452832  16
length(detP[detP > 0.01])  # 3808
max(colMeans(detP))  # 0.0001333676
detP[1:3, 1:5]
#            S1.baseline S1.loading S2.baseline S2.loading S3.baseline
# cg00000029           0          0           0          0           0
# cg00000109           0          0           0          0           0
# cg00000165           0          0           0          0           0


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
#  1340 451492

keep <- rowSums(detP < 0.01) > 0.95*ncol(mlset)
mlset <- mlset[keep,]; dim(mlset)  # 451492  16
detP <- detP[keep,]; dim(detP)  # 451492  16



#### removing outlier (b & M vals) ####
mVals.pre <- estimateM(as(mlset, "MethyLumiM"), returnType = "matrix")
bVals.pre <- mlset@assayData$betas

par(mfrow=c(1,1))
det.outlier.m.pre <- detectOutlier(mVals.pre, ifPlot = TRUE)  # S1.baseline
det.outlier.b.pre <- detectOutlier(bVals.pre, ifPlot = TRUE)  # S1.baseline

#### mlset からoutlierをremove
mlset <- mlset[,colnames(mlset) != "S1.baseline"]; dim (mlset)  # 451492       15
pD <- pD[rownames(pD) != "S1.baseline",]



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

par(mfrow=c(1,3))
plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$timepoint)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$sample)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$sample)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$slide)], cex=0.8)
legend("bottomright", legend=levels(factor(pD$slide)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")


#. Examine Higher dimensions ----
par(mfrow=c(1,3))
plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], dim=c(1,3))
legend("bottomright", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], dim=c(2,3))
legend("bottomleft", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")

plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], dim=c(3,4))
legend("topleft", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")


#. density plot: M & beta values ----
par(mfrow=c(1,2))
plotDensity(mVals, legend = FALSE, main = "M values", col=pal[factor(pD$timepoint)])
legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")
plotDensity(bVals, legend = FALSE, main = "Beta values", col=pal[factor(pD$timepoint)])
legend("topleft", legend=levels(factor(pD$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")



#### ComBat: removing batch effect ####
# mod <- model.matrix(~ as.factor(pD$timepoint) + as.factor(pD$sample), data = pD)
# mod0 <- model.matrix(~ as.factor(pD$sample), data = pD)

batch <- as.factor(pD$slide)
modcombat <- model.matrix(~1, data = pD)
combat.edata <- ComBat(dat = mVals, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = TRUE)
# pValuesComBat <- f.pvalue(combat.edata, mod, mod0)
# qValuesComBat <- p.adjust(pValuesComBat, method="BH")


#. MDS: visualizing Combat()ed data ----
par(mfrow=c(1,3))
plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pD$timepoint)], cex=0.8)
legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pD$sample)], cex=0.8)
legend("topright", legend=levels(factor(pD$sample)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pD$slide)], cex=0.8)
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
contMatrix <- makeContrasts(loading - baseline, levels=design)
# contMatrix

#. Fit the Contrasts ----
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

#. Numbers of DM CpGs at FDR < 0.05 ----
summary(decideTests(fit2))
summary(decideTests(fit2, p.value = 0.2))
summary(decideTests(fit2, adjust.method = "none", p.value = 0.05)) # p-value
# loading - baseline (FDR < 0.05)
# Down                    0
# NotSig             451492
# Up                      0

# loading - baseline (FDR < 0.2)
# Down                    0
# NotSig             451492
# Up                      0

# loading - baseline (p < 0.05)
# Down                16678
# NotSig             423662
# Up                  11152




#### Get the Table of Results ####
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(combat.edata), ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
# anno <- AnnotationDbi::select(hgu133plus2.db,
#                               keys = (ann450kSub$UCSC_RefGene_Name),
#                               columns = c("ENSEMBL", "ENTREZID"),
#                               keytype = "PROBEID")
# head(ann450kSub)

DMPs <- topTable(fit2, num=Inf, coef=1, genelist = ann450kSub, sort.by = "p")
DMPs[1:3,c("logFC", "P.Value", "adj.P.Val")]
#                 logFC      P.Value adj.P.Val
# cg01500115  0.5109598 6.032010e-06 0.7093039
# cg06081361 -0.4651828 1.135276e-05 0.7093039
# cg03969906  0.5554859 1.169004e-05 0.7093039

DMPs %>% select(logFC, P.Value, adj.P.Val) %>% 
    filter(P.Value < 0.05 & abs(logFC) > 1) %>% 
    nrow()  # 3

summary(DMPs$logFC)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -1.704712 -0.072202 -0.010060 -0.009545  0.052321  1.848543


write.table(DMPs, file="DMPs.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0905.csv", sep=",", row.names=FALSE)
dmps.005 <- subset(DMPs, P.Value < 0.05)
write.table(dmps.005, file="DMPs.p005.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv", sep=",", row.names=FALSE)



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

labs(title = "GSE114763 Resistance p < 0.05", titleLabSize = 10,
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
cpgs.sig.names <- dmps.005$Name; length(cpgs.sig.names)   # 27830
cpgs.pos.names <- dmps.005[dmps.005$logFC > 0,]$Name; length(cpgs.pos.names)  # 11152
cpgs.neg.names <- dmps.005[dmps.005$logFC < 0,]$Name; length(cpgs.neg.names)  # 16678
all.names <- DMPs$Name; length(all.names)  # 451492
gst <- gometh(sig.cpg=cpgs.sig.names, all.cpg=all.names, plot.bias=F)
gst.pos <- gometh(sig.cpg=cpgs.pos.names, all.cpg=all.names, plot.bias=F)
gst.neg <- gometh(sig.cpg=cpgs.neg.names, all.cpg=all.names, plot.bias=F)

#. categories ----
t10GO <- topGSA(gst, num=Inf)
t10GO <- t10GO[t10GO$P.DE < 0.05,]; length(t10GO$TERM)  # 902
t10GO.pos <- topGSA(gst.pos, num=Inf)
t10GO.pos <- t10GO.pos[t10GO.pos$P.DE < 0.05,]; length(t10GO.pos$TERM)  # 464
t10GO.neg <- topGSA(gst.neg, num=Inf)
t10GO.neg <- t10GO.neg[t10GO.neg$P.DE < 0.05,]; length(t10GO.neg$TERM)  # 1113

t10GO$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO$TERM)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO$TERM)  # 0
t10GO.pos$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)]
grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)  # 0
t10GO.neg$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.neg$TERM)]  # 0
# [1] "regulation of blood pressure"
grep("blood pressure|cardio vascular|inflamation", t10GO.neg$TERM)  # 847


write.csv(t10GO, "go.p005.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")
write.csv(t10GO.pos, "go.p005.pos.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")
write.csv(t10GO.neg, "go.p005.neg.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")
t10GO$TERM[1:20]
# [1] "anatomical structure morphogenesis"             "cell adhesion"                                 
# [3] "biological adhesion"                            "cell junction"                                 
# [5] "anchoring junction"                             "positive regulation of myotube differentiation"
# [7] "embryonic skeletal system morphogenesis"        "nervous system development"                    
# [9] "movement of cell or subcellular component"      "cell-cell adhesion"                            
# [11] "muscle contraction"                             "muscle system process"                         
# [13] "myofibril"                                      "muscle structure development"                  
# [15] "animal organ morphogenesis"                     "cellular component morphogenesis"              
# [17] "cell surface receptor signaling pathway"        "peptidyl-tyrosine modification"                
# [19] "contractile fiber"                              "regulation of developmental process"

