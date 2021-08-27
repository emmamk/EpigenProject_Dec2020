#### 450k DMA ####
# Endurance exercise: Before, After
# 参照URL: https://kasperdanielhansen.github.io/genbioconductor/html/minfi.html
# set up a path to the data directory
# setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/GSE60655_End')
dataDir <- ('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/DataDir/')

##### Obtaining the Data #####
library(GEOquery)
getGEOSuppFiles("GSE60655", makeDirectory = F) # 450k
untar("GSE60655_RAW.tar", exdir = "idat")
head(list.files("idat", pattern = "idat"))

# idatFiles <- list.files("idat", pattern = "idat.gz$", full = TRUE)
idatFiles <- list.files("idat", pattern = "idat", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

##### Libraries #####
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(dplyr)
#library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
#library(IlluminaHumanMethylationEPICmanifest)
#library(BSgenome.Hsapiens.UCSC.hg19)

##### Get 450k Annotation Data #####
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
head(ann450k)

##### Read Raw Data From IDAT Files ####
# (Ignore Warnings) 
rgSet <- read.metharray.exp("idat")
rgSet

pData(rgSet)
head(sampleNames(rgSet))

##### Rename Samples with Descriptive Ones #####
# library(GEOquery)
geoMat <- getGEO("GSE60655", getGPL = FALSE)
pD.all <- pData(geoMat[[1]])
head(pD.all)

# pD.all$treatment_protocol_ch1
# Healthy, sedentary subjects performed supervised one-legged, knee-extension exercise training for three months.
# Skeletal muscle biopsies from the vastus lateralis were taken at rest, before and after the training period.

# pD.all$data_processing
# Raw intesities exported from GenomeStudio were imported into R.
# Probes on chromosomes X and Y were removed.
# Probes with a detection P > 0.01 exceeding 5% of the samples were also filtered out. 
# Color-bias adjustment and quantile normalization were performed as implemented in lumi.
# BMIQ correction was applied to correct for probe type bias.
# Finally, M-values were batch-corrected usign the ComBat function.

pD <- as(pD.all[, c("title", "gender:ch1", "source_name_ch1", "subject:ch1", "geo_accession")], "DataFrame")
# (as( ,"DataFrame"))が無いと後の pData(rgSet) <- pD がerror.
# pData(rgSet)はDataFrameである必要がある、とdeveloperが修正している:
# https://github.com/genomicsclass/labs/pull/90
# 参考: https://github.com/hansenlab/minfi/issues/174
head(pD)
head(pD.all)

# name: don't start with numbers nor contain "-"
names(pD)[c(1,2,3,4)] <- c("rep", "gender", "timepoint", "sample")
pD$rep <- str_extract(pD$rep, pattern="rep.")
pD$timepoint <- str_extract(pD$timepoint, pattern="before|after")
pD$sample <- paste0("S", pD$sample)
pD$ID <- paste(pD$sample, pD$timepoint, sep=".")
head(pD)
tail(pD)

#### Tidy Up the Data ####
#. Remove Replicates ----
keep <- grep("rep1", pD$rep)
pD <- pD[keep,]
rgSet <- rgSet[,keep]

#. Create Subset ----
keep <- grep("^male$", pD$gender)
pD <- pD[keep,]
rgSet <- rgSet[,keep]

dim(pD) # 14  6
dim(rgSet)  # 622399     14
table(pD$gender)
# male 
# 14 


# sample:S8はReloadのデータが無いのでdel.
# pD <- pD[-15,]
# rgSet <- rgSet[,-15]

#. Merge Target Data & rgSet ----
sampleNames(rgSet) <- str_sub(sampleNames(rgSet), 1, 10)
sampleNames(rgSet) <- rownames(pD) <- pD$ID[match(sampleNames(rgSet), pD$geo_accession)]
pD <- pD[sampleNames(rgSet),]
pData(rgSet) <- pD
rgSet

#### Quality Control ####
#. Detection p-values ----
detP <- detectionP(rgSet)
detP[1:3, 1:5]
# S1.before      S1.after    S4.before      S4.after    S5.before
# cg00050873  0.000000e+00  0.000000e+00  0.00000e+00  0.000000e+00 0.000000e+00
# cg00212031  0.000000e+00  0.000000e+00  0.00000e+00  0.000000e+00 0.000000e+00
# cg00213748 3.333207e-103 2.922029e-120 3.99785e-104 3.585549e-146 2.210492e-85
max(colMeans(detP))  # 0.0003285641


#. Plot detP ----
# examine mean detection p-values across all samples to identify any failed samples
pal <- brewer.pal(8,"Dark2")
# par(mfrow=c(1,1))

# All・Timepoint
barplot(colMeans(detP), col=pal[factor(pD$timepoint)], las=2, 
        cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")
# legend("topleft", legend=levels(factor(pD$timepoint)), fill=pal,
#       bg="white", cex = 0.5)

#. qcReport(minfi) ----
# Ex.of another useful quality control plot
# 450K用? 850Kは不明
qcReport(rgSet, sampNames=pD$ID, sampGroups=pD$timepoint, 
         pdf="qcReport_male.pdf")

#. Remove Poor Quality Samples ----
max(colMeans(detP))
keep <- colMeans(detP) < 0.05
rgSet <- rgSet[,keep]
rgSet

# remove from "target" data
pD <- pD[keep,]
pD[,1:5]

# remove from detection p-value table
detP <- detP[,keep]
dim(detP)  # 485512     14



#### Normalisation ####
# normalize the data; this results in a GenomicRatioSet object
# choose optimal normalization method suitable for the data
mSetSq <- preprocessQuantile(rgSet)

# create a MethylSet object from the raw data for later plotting
mSetRaw <- preprocessRaw(rgSet)

#. Visualise Before & After Normalisation ----
par(mfrow=c(1,2))
# Before
densityPlot(rgSet, sampGroups=pD$timepoint, main="Raw", legend=FALSE)
# legend("top", legend = levels(factor(pD$timepoint)),
#       text.col=brewer.pal(8,"Dark2"), cex=0.5)

# After
densityPlot(getBeta(mSetSq), sampGroups=pD$timepoint, main="Normalized", legend=FALSE)
# legend("top", legend = levels(factor(pD$timepoint)),
#       text.col=brewer.pal(8,"Dark2"), cex=0.5)

#. Data exploration ----
#. . MDS plot ----
# Multi-dimensional scaling (MDS) plots: to look at largest sources of variation
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], cex=0.6)
# legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal,
#      fill=pal, bg="white", cex=0.5)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",  
        col=pal[factor(pD$sample)], cex=0.6)
# legend("topright", legend=levels(factor(pD$sample)), text.col=pal,
#      fill=pal, bg="white", cex=0.5)

# plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
# col=pal[factor(pD$gender)], cex=0.6)

#. . Examine Higher dimensions ----
par(mfrow=c(1,3))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], dim=c(1,3))
# legend("top", legend=levels(factor(pD$timepoint)), text.col=pal,
#       cex=0.5, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], dim=c(2,3))
# legend("topleft", legend=levels(factor(pD$timepoint)), text.col=pal,
#       cex=0.5, bg="white")

plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], dim=c(3,4))
# legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal,
#       cex=0.5, bg="white")

#. Filtering ----
#. . Poor Performing Probes ----
# Poor performing probes(unreliable signal) are generally filtered out
# prior to differential methylation analysis.
# i.e. fewer statistical tests

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

# remove any probes that have failed in one or more samples
# ex.detPのcolに1つNAがあるとrowSums(detP)は9、一方ncol(mSetSq)は常に10.
# そのように左右の数が合わない物を除外する.
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
dim(mSetSqFlt)  # 484222     14


#. . Probes on the Sex Chromosomes ----
# Only applies if the dataset has both gender
# keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# table(keep)
# mSetSqFlt <- mSetSqFlt[keep,]

#. . Probes with SNPs at CpG Site ----
# (common SNPs may affect the CpG)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
dim(mSetSqFlt)  # 467244     14


#. . Cross Reactive Probes ----
# (multiple places in the genome)
# Chen et al. 2013: 450k
# https://www.tandfonline.com/doi/full/10.4161/epi.23470

#. . . Obtaining Reference Files for 450k ----
# https://bioc.ism.ac.jp/packages/devel/bioc/vignettes/ramwas/inst/doc/RW5a_matrix.html#probes-with-snps-and-in-cross-reactive-regions
host = "https://wapps.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/"
files = c(
    i450k_ns.xlsx = "48639-non-specific-probes-Illumina450k.xlsx",
    i450k_pl.xlsx = "48640-polymorphic-CpGs-Illumina450k.xlsx")

for( i in seq_along(files) )
    download.file(
        url = paste0(host, files[i]),
        destfile = names(files)[i],
        mode = "wb",
        quiet = TRUE)

#. . . Ref Files ----
library(readxl)
ex1 = read_excel(paste(dataDir, "i450k_ns.xlsx", sep="/"), sheet = 1)
ex2 = read_excel(paste(dataDir, "i450k_pl.xlsx", sep="/"), sheet = 1)
ex3 = read_excel(paste(dataDir, "i450k_pl.xlsx", sep="/"), sheet = 2)

exclude.snp = unique(c(
    ex1$TargetID,
    ex2$PROBE,
    ex3$PROBE[ (ex3$BASE_FROM_SBE < 10) & (ex3$AF > 0.01)]))

keep <- !(featureNames(mSetSqFlt) %in% exclude.snp)
table(keep)
# keep
# FALSE   TRUE 
# 82024 385220
mSetSqFlt <- mSetSqFlt[keep,] 
dim(mSetSqFlt)  # 385220     14

# rm(host, files, i, ex1, ex2, ex3)
rm(ex1, ex2, ex3)

#. Re-examine MDS Plots ----
# Once the data has been filtered and normalised, it is often useful to re-examine
# the MDS plots to see if the relationship between the samples has changed
par(mfrow=c(1,1))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], cex=0.65)
# legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal,
#       cex=0.5, bg="white")

#. . Before & After Removing Gender-probes ----
par(mfrow=c(1,2))
# Before
plotMDS(getM(mSetSq), top=1000, gene.selection="common", main="Before", 
        col=pal[factor(pD$gender)], cex=0.6)

# After
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", main="After",
        col=pal[factor(pD$gender)], cex=0.65)

#. . Remove Female Samples ----
# keep <- pD$gender=="male"
# table(keep)
# pD <- pD[keep,]
# mSetSqFlt <- mSetSqFlt[,keep]
# mSetSqFlt

#### Calculate M & beta-values ####
# M-values(getM(minfi)):
# nicer statistical properties: better for use in statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])
# S1.before  S1.after S4.before  S4.after S5.before
# cg13869341  2.505376  2.566522  2.515211  2.719321  2.359755
# cg24669183  2.182828  2.000885  1.839788  2.179010  1.814513
# cg15560884  1.384137  1.322324  1.414610  1.109830  1.292936
# cg01014490 -3.950538 -4.716351 -5.211276 -3.837747 -4.396400
# cg17505339  3.156644  3.528931  3.490768  3.661623  3.925537
# cg11954957  2.154270  1.776505  1.862762  2.067354  1.935827


# beta-values(getbeta(minfi)):
# are easy to interpret: better for displaying data
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])
# S1.before  S1.after  S4.before   S4.after  S5.before
# cg13869341 0.8502540 0.8555705 0.85111982 0.86817177 0.83694322
# cg24669183 0.8195085 0.8000982 0.78163912 0.81911669 0.77863420
# cg15560884 0.7230038 0.7143417 0.72721395 0.68336408 0.71016695
# cg01014490 0.0607506 0.0366457 0.02628343 0.06536781 0.04533192
# cg17505339 0.8991684 0.9202735 0.91831096 0.92676616 0.93825303
# cg11954957 0.8165620 0.7740599 0.78434483 0.80736534 0.79278794


#. Visualise M & beta-values ----
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=pD$timepoint, main="Beta values", 
            legend=FALSE, xlab="Beta values")
# legend("top", legend = levels(factor(pD$timepoint)),
#       text.col=brewer.pal(8,"Dark2"), cex=0.5)

densityPlot(mVals, sampGroups=pD$timepoint, main="M-values", 
            legend=FALSE, xlab="M values")
# legend("topleft", legend = levels(factor(pD$timepoint)),
#       text.col=brewer.pal(8,"Dark2"), cex=0.5)

#### Probe-wise Differential Methylation Analysis ####
#. Create a Design Matrix ----
# factor of interest
timepoint <- factor(pD$timepoint)
relevel(timepoint, "before")
# individual effect that we need to account for
individual <- factor(pD$sample)

# design matrix
design <- model.matrix(~0+timepoint+individual, data=pD)
colnames(design) <- c(levels(timepoint),levels(individual)[-1])
design[1:5, 1:5]
# (~の後にグループ名を入れる. 0が無いとグループ名の最後がrefになる.
# relevelでrefの変更が可能.
# 参照: https://genomicsclass.github.io/book/pages/expressing_design_formula.html

#. Fit the Linear Model ----
# lmFit: Linear Model For Series Of Arrays
fit <- lmFit(mVals, design)
head(mVals)
head(fit)

#. Contrast Matrix ----
# for specific comparisons
contMatrix <- makeContrasts(before-after, levels=design)
contMatrix

#. Fit the Contrasts ----
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)
head(fit2)

#. Numbers of DM CpGs at FDR < 0.05 ----
summary(decideTests(fit2))
summary(decideTests(fit2, p.value = 0.2))
summary(decideTests(fit2, adjust.method = "none", p.value = 0.05)) # p-value
# before - after
# Down            16459
# NotSig         352598
# Up              16163
summary(decideTests(fit2, adjust.method = "none", p.value = 0.01))
# before - after
# Down             2911
# NotSig         379033
# Up               3276


#. Get the Table of Results ----
# for the first contrast: coef = 1
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
head(ann450kSub)

DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)
# (DMPs = differentially methylated probes )
# (To order by p-value, the user can specify topTable(sort.by="p").
# DefultはB-statistic("B"). 但しmost casesでp値とB-statisticはidentical.)
# num=Inf はデータを全て摘出. Cut Offを設けていない.

##### Explore Results #####
# summary(DMPs)
# table(DMPs$adj.P.Val < 0.2)  # 0
# sum(DMPs$adj.P.Val < 0.1, na.rm=TRUE)
# resSig_P005 <- subset(DMPs, P.Value < 0.05)
# nrow(resSig_P005)
# head(resSig_P005[ order( resSig_P005$logFC ), ])
# head(resSig_P005[ order( -resSig_P005$logFC ), ])
# resOrdered_P005 <- resSig_P005[order(resSig_P005$P.Value),]
# head(select(resOrdered_P005, P.Value))

write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)
write.table(resOrdered_P005, file="DMPs_p005.csv", sep=",", row.names=FALSE)
# 数百MB, 開けるのに時間かかる.

#. Plot Top 4 Most Sig DM CpGs ----
# 出力したいデータ = 第一引数
par(mfrow=c(2,2))
sapply(rownames(resOrdered_P005)[1:4], function(cpg){
    plotCpg(bVals, cpg=cpg, pheno=pD$timepoint, ylab = "Beta values")
})

#### Differential Methylation Analysis of Regions ####
# bumphunter(minfi), dmrFind(charm): slow
# dmrcate(DMRcate): faster, based on limma = can use design & contMatrix assigned earlier
# memo. cpg.annotate(): Annotate CpGs With Their Chromosome Position And Test Statistic
myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "before - after",
                             arraytype = "450K")
str(myAnnotation)

# dmrcate(): identify significantly differentially (or variable) methylated regions.
# ?dmrcate(): # pcutoffのデフォルトは"fdr"
DMRs <- dmrcate(myAnnotation, pcutoff = 0.05, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges

#. Visualise the Results ----
# to ensure that they make sense.
# set up the grouping variables and colours
groups <- pal[1:length(unique(pD$timepoint))]
names(groups) <- levels(factor(pD$timepoint))
cols <- groups[as.character(factor(pD$timepoint))]

# draw the plot for the top DMR
par(mfrow=c(1,1))
DMR.plot(ranges = results.ranges, dmr = 1, CpGs = bVals, phen.col = cols, 
         what = "Beta", arraytype = "450K", genome = "hg19")

#. Customising visualisations of methylation data ----
### 122120: やってない ###
# indicate which genome is being used
gen <- "hg19"
# the index of the DMR that we will plot 
dmrIndex <- 1
# extract chromosome number and location from DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))
# add 25% extra space to plot
minbase <- start - (0.25*(end-start))
maxbase <- end + (0.25*(end-start))

# CpG islands
# cf. "/model-based-cpg-islands-hg19-chr17.txt"はlibrary(BSgenome.Hsapiens.UCSC.hg19)と同一
# 拡大して実際のgene配列を見ることもできる
# 参考: http://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html#3_Plotting_parameters
islandHMM <- read.csv(paste0(dataDirectory,
                             "/model-based-cpg-islands-hg19-chr17.txt"),
                      sep="\t", stringsAsFactors=FALSE, header=FALSE)
head(islandHMM)
# (V1:V8の各column names: chr, start, end, length, CpGcount, GCcontent, pctGC, obsExp)
# http://www.haowulab.org/software/makeCGI/model-based-cpg-islands-hg19.txt

# Rle(): Run Length Encoding
# strand(*): the exact strand of the location is unknown
islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
islandData

# DNAseI hypersensitive sites
# wgEncodeRegDnaseClusteredV3chr17.bed: UCSCの下記ブラウザーでDL可能
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=977545497_nv3CH2qiBjrVdLZ82EACf7OLzjwI&clade=mammal&org=Human&db=hg19&hgta_group=regulation&hgta_track=wgEncodeRegDnaseClustered&hgta_table=0&hgta_regionType=range&position=chrX%3A15%2C578%2C261-15%2C621%2C068&hgta_outputType=primaryTable&hgta_outFileName=
dnase <- read.csv(paste0(dataDirectory,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                  sep="\t",stringsAsFactors=FALSE,header=FALSE)
head(dnase)

dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))),
                     data=dnase[,5])
dnaseData

# set up the ideogram, genome and RefSeq tracks
# that will provide context for our methylation data.
iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue",stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE)

# Ensure that the methylation data is ordered by
# chromosome and base position.
annEPICOrd <- annEPICSub[order(annEPICSub$chr,annEPICSub$pos),]
head(annEPICOrd)

bValsOrd <- bVals[match(annEPICOrd$Name,rownames(bVals)),]
head(bValsOrd)

# create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(annEPICOrd$chr),
                   ranges=IRanges(start=annEPICOrd$pos, end=annEPICOrd$pos),
                   strand=Rle(rep("*",nrow(annEPICOrd))),
                   betas=bValsOrd)
# extract data on CpGs in DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])

# methylation data track
methTrack <- DataTrack(range=cpgData, groups=pD$timepoint,genome = gen,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)
# CpG island track
islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                               chromosome=chrom,fill="darkgreen")

# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)

# DMR position data track
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom, fill="darkred")

# Set up the track list and indicate
# the relative sizes of the different tracks.
tracks <- list(iTrack, gTrack, methTrack, dmrTrack,
               islandTrack, dnaseTrack, rTrack)
sizes <- c(2,2,5,2,2,2,3)

# Draw the plot using the plotTracks function.
plotTracks(tracks, from=minbase, to=maxbase,
           showTitle=TRUE, add53=TRUE, add35=TRUE,
           grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
# (lty = line type: 0. “blank”, 1. “solid”, 2. “dashed”, 3. “dotted”,
#                   4. “dotdash”, 5. “longdash” and 6. “twodash”.)

#### Additional analyses ####
#. Gene Ontology Testing ####
# Get the significant CpG sites at less than 5% FDR
# 010321: GSE114763とのOverlapsを抽出
head(DMPs)
GSE60655Sig <- subset(DMPs, P.Value < 0.05)
GSE60655_P <- subset(GSE60655Sig, logFC > 0)
GSE60655_N <- subset(GSE60655Sig, logFC < 0)

GSE114763Sig <- read.csv("GSE114763/DMPs.csv")
GSE114763Sig <- subset(GSE114763Sig, P.Value < 0.05)
GSE114763_P <- subset(GSE114763Sig, logFC > 0)
GSE114763_N <- subset(GSE114763Sig, logFC < 0)

# intersect()よりもinner_join()で他の変数も残した方が良さそう
# hyperとhypoどちらも一緒にGO解析に渡して良いか検討中.
# OLCpGs_P <- intersect(GSE60655_P$Name, GSE114763_P$Name)
# OLCpGs_N <- intersect(GSE60655_N$Name, GSE114763_N$Name)

# length(OLCpGs_P)
# length(OLCpGs_N)

# OverlapCpGsを解析するためにマスターコードを変更(12/26/20)
# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
length(all)

par(mfrow=c(1,2))
gst_P <- gometh(sig.cpg=OLCpGs_P, all.cpg=all, plot.bias=TRUE)
gst_N <- gometh(sig.cpg=OLCpGs_N, all.cpg=all, plot.bias=TRUE)
# warning regarding multiple symbols will always be displayed as
# there are genes that have more than one alias; not a cause for concern.

#. . Top 10 GO Categories ----
t10GO_P <- topGSA(gst_P, num=Inf)
t10GO_P <- t10GO_P[t10GO_P$P.DE < 0.05,]
t10GO_N <- topGSA(gst_N, num=Inf)
t10GO_N <- t10GO_N[t10GO_N$P.DE < 0.05,]
length(t10GO_P$TERM)
length(t10GO_N$TERM)

t10GO_P$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO_P$TERM)]
t10GO_N$TERM[grep("blood pressure|cardio vascular", t10GO_N$TERM)]

t10GO <- rbind(t10GO_P, t10GO_N)
write.csv(t10GO, "go.csv")
