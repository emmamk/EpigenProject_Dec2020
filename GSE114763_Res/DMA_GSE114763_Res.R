#### EPIC DMA ####
# Resistance exercise: Base, Right After, 7wLoading, 7wUnloading, Reloading
# 参照URL: https://kasperdanielhansen.github.io/genbioconductor/html/minfi.html
# setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/')
dataDir <- ('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/DataDir/')

##### Obtaining the Data #####
library(GEOquery)
getGEOSuppFiles("GSE114763", makeDirectory = F) # EPIC
untar("GSE114763_RAW.tar", exdir = "idat")
head(list.files("idat", pattern = "idat"))

idatFiles <- list.files("idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)

##### Libraries #####
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(dplyr)
#library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#library(IlluminaHumanMethylation450kmanifest)
#library(BSgenome.Hsapiens.UCSC.hg19)

##### Get EPIC Annotation Data #####
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
head(annEPIC)

##### Read Raw Data From IDAT Files ####
# (Ignore Warnings) 
rgSet <- read.metharray.exp("idat")
rgSet

pData(rgSet)
head(sampleNames(rgSet))

##### Rename Samples with Descriptive Ones #####
# library(GEOquery)
geoMat <- getGEO("GSE114763")
pD.all <- pData(geoMat[[1]])
head(pD.all)

pD <- as(pD.all[, c("title", "title", "geo_accession")], "DataFrame")
# (as( ,"DataFrame"))が無いと後の pData(rgSet) <- pD がerror.
# pData(rgSet)はDataFrameである必要がある、とdevelopperが修正している:
# https://github.com/genomicsclass/labs/pull/90
# 参考: https://github.com/hansenlab/minfi/issues/174
head(pD)
head(pD.all)

# name: don't start with numbers nor contain "-"
names(pD)[c(1,2)] <- c("individual", "timepoint")
pD$individual <- str_extract(pD$individual, pattern="._Rep.")
pD$timepoint <- gsub("SkM_Epi_Mem_.*: Muscle_|7wk_|_P.*", "", pD$timepoint)
pD$timepoint <- str_sub(pD$timepoint, 1, 2)
pD$sample <- paste0("S", str_sub(pD$individual, 1, 1))
pD$ID <- paste(pD$timepoint, str_sub(pD$individual, 1, 1), sep=".")
head(pD)
tail(pD)

#### Tidy Up the Data ####
#. Remove Replicates ----
keep <- grep("Rep1", pD$individual)
pD <- pD[keep,]
rgSet <- rgSet[,keep]

#. Create Subset ----
# keep <- grep("Ba|Re", pD$timepoint)
# pD <- pD[keep,]
# rgSet <- rgSet[,keep]

# sample:S8はReloadのデータが無いのでdel.
pD <- pD[-15,]
rgSet <- rgSet[,-15]

#. Merge Target Data & rgSet ----
sampleNames(rgSet) <- rownames(pD) <- pD$ID
pD <- pD[sampleNames(rgSet),]
pData(rgSet) <- pD
rgSet

#### Quality Control ####
#. Detection p-values ----
detP <- detectionP(rgSet)
max(colMeans(detP))
head(detP)

#. Plot detP ----
pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,1))

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
         pdf="qcReport.pdf")

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
dim(detP)

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
#        text.col=brewer.pal(8,"Dark2"), cex=0.5)

# After
densityPlot(getBeta(mSetSq), sampGroups=pD$timepoint,
            main="Normalized", legend=FALSE)
# legend("top", legend = levels(factor(pD$timepoint)), 
#        text.col=brewer.pal(8,"Dark2"), cex=0.5)

#. Data exploration ----
#. . MDS plot ----
# Multi-dimensional scaling (MDS) plots: to look at largest sources of variation
par(mfrow=c(1,1))
plotMDS(getM(mSetSq), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], cex=0.6)
# legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal,
#      fill=pal, bg="white", cex=0.5)

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
mSetSqFlt

#. . Probes on the Sex Chromosomes ----
# Only applies if the dataset has both gender
# keep <- !(featureNames(mSetSqFlt) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
# table(keep)
# mSetSqFlt <- mSetSqFlt[keep,]

#. . Probes with SNPs at CpG Site ----
# (common SNPs may affect the CpG)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

#. . Cross Reactive Probes ----
# (multiple places in the genome)
# Pidsley et al. 2016: EPIC
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1#Sec22

#. . . Obtaining Reference Files for EPIC ----
# https://bioc.ism.ac.jp/packages/devel/bioc/vignettes/ramwas/inst/doc/RW5a_matrix.html#probes-with-snps-and-in-cross-reactive-regions
host = "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/"
files = c(
    S1_cross_reactive.csv     = '13059_2016_1066_MOESM1_ESM.csv',
    S4_snp_cpg.csv            = '13059_2016_1066_MOESM4_ESM.csv',
    S5_snp_base_extension.csv = '13059_2016_1066_MOESM5_ESM.csv',
    S6_snp_body.csv           = '13059_2016_1066_MOESM6_ESM.csv')

for( i in seq_along(files) )
    download.file(
        url = paste0(host, files[i]),
        destfile = names(files)[i],
        mode = "wb",
        quiet = TRUE)

#. . . Ref Files ----
snpcpgs1 = read.csv(paste(dataDir,'S1_cross_reactive.csv', sep="/"), stringsAsFactors = FALSE)
snpcpgs4 = read.csv(paste(dataDir,'S4_snp_cpg.csv', sep="/"), stringsAsFactors = FALSE)
snpcpgs5 = read.csv(paste(dataDir,'S5_snp_base_extension.csv', sep="/"), stringsAsFactors = FALSE)
snpcpgs6 = read.csv(paste(dataDir,'S6_snp_body.csv', sep="/"), stringsAsFactors = FALSE)

xReactiveProbes = unique(c(
    snpcpgs1$X,
    snpcpgs4$PROBE,
    snpcpgs5$PROBE,
    snpcpgs6$PROBE[
        pmax(snpcpgs6$VARIANT_START - snpcpgs6$MAPINFO, 
             snpcpgs6$MAPINFO - snpcpgs6$VARIANT_END) < 10]))

keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt

rm(snpcpgs1, snpcpgs4, snpcpgs5, snpcpgs6)

#. Re-examine MDS Plots ----
# Once the data has been filtered and normalised, it is often useful to re-examine
# the MDS plots to see if the relationship between the samples has changed
# ?plotMDS()
par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], cex=0.65)
legend("bottomleft", legend=levels(factor(pD$timepoint)), text.col=pal,
      cex=0.5, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(pD$sample)], cex=0.65)
legend("bottomleft", legend=levels(factor(pD$timepoint)), text.col=pal,
       cex=0.5, bg="white")
#. . Examine Higher dimensions ----
par(mfrow=c(1,3))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], dim=c(1,3))
# legend("top", legend=levels(factor(pD$timepoint)), text.col=pal,
#       cex=0.5, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], dim=c(2,3))
# legend("topleft", legend=levels(factor(pD$timepoint)), text.col=pal,
#       cex=0.5, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
        col=pal[factor(pD$timepoint)], dim=c(3,4))
# legend("topright", legend=levels(factor(pD$timepoint)), text.col=pal,
#       cex=0.5, bg="white")

# ここで良いのか考え中 ----
# プロットで見た結果、Un.8が飛び抜けてはみ出しているので削除
pD <- pD[-39,]
rgSet <- rgSet[,-39]
mSetSqFlt <- mSetSqFlt[,-39] 
mSetSqFlt

#### Calculate M & beta-values ####
# M-values(getM(minfi)):
# nicer statistical properties: better for use in statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

# beta-values(getbeta(minfi)):
# are easy to interpret: better for displaying data
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

#. Visualise M & beta-values ----
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=pD$timepoint, main="Beta values", 
            legend=FALSE, xlab="Beta values")
legend("top", legend = levels(factor(pD$timepoint)), 
       text.col=brewer.pal(8,"Dark2"), cex=0.5)

densityPlot(mVals, sampGroups=pD$timepoint, main="M-values", 
            legend=FALSE, xlab="M values")
legend("topleft", legend = levels(factor(pD$timepoint)), 
       text.col=brewer.pal(8,"Dark2"), cex=0.5)

#### Probe-wise Differential Methylation Analysis ####
#. Create a Design Matrix ----
# factor of interest
Timepoint <- factor(pD$timepoint)
Timepoint <- relevel(Timepoint, "Ba", "Ac", "Lo", "Un", "Re")
# individual effect that we need to account for
Individual <- factor(pD$sample)

# design matrix
# ?model.matrix()
design <- model.matrix(~0+Timepoint+Individual, data=pD)
colnames(design)
colnames(design) <- c(levels(Timepoint),levels(Individual)[-1])
# (~の後にグループ名を入れる. 0が無いとグループ名の最後がrefになる.
# relevelでrefの変更が可能.
# 参照: https://genomicsclass.github.io/book/pages/expressing_design_formula.html

#. Fit the Linear Model ----
# lmFit: Linear Model For Series Of Arrays
# ?lmFit()
mVals <- mVals
fit <- lmFit(mVals, design)

#. Contrast Matrix ----
# for specific comparisons
contMatrix <- makeContrasts(Ba-Un,
                            Un-Re,
                            Ba-Ac,
                            Ba-Lo,
                            Ba-Re,
                            levels=design)
contMatrix

#. Fit the Contrasts ----
# ?contrasts.fit()
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

#. Numbers of DM CpGs at FDR < 0.05 ----
summary(decideTests(fit2))

#. Get the Table of Results ----
# for the first contrast: coef = 1
annEPICSub <- annEPIC[match(rownames(mVals),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))]

DMPs_BaUn <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
DMPs_UnRe <- topTable(fit2, num=Inf, coef=2, genelist=annEPICSub)
head(DMPs)
# (DMPs = differentially methylated probes )
# (To order by p-value, the user can specify topTable(sort.by="p").
# DefultはB-statistic("B"). 但しmost casesでp値とB-statisticはidentical.)
# num=Inf はデータを全て摘出. Cut Offを設けていない.

##### Explore Results #####
summary(DMPs_BaUn)
summary(DMPs_UnRe)
table(DMPs_BaUn$adj.P.Val < 0.1)
table(DMPs_UnRe$adj.P.Val < 0.1)
sum(DMPs_BaUn$adj.P.Val < 0.1, na.rm=TRUE)
sum(DMPs_UnRe$adj.P.Val < 0.1, na.rm=TRUE)
resSig_P005BaUn <- subset(DMPs_BaUn, P.Value < 0.05)
resSig_P005UnRe <- subset(DMPs_UnRe, P.Value < 0.05)
dim(resSig_P005BaUn)
dim(resSig_P005UnRe)
head(resSig_P005[ order( resSig_P005$logFC ), ])
head(resSig_P005[ order( -resSig_P005$logFC ), ])
resOrdered_P005BaUn <- resSig_P005BaUn[order(resSig_P005BaUn$P.Value),]
resOrdered_P005UnRe <- resSig_P005UnRe[order(resSig_P005UnRe$P.Value),]
head(select(resOrdered_P005, P.Value))

write.table(resOrdered_P005BaUn, file="DMPs_Res_P005BaUn.csv", sep=",", row.names=FALSE)
write.table(resOrdered_P005UnRe, file="DMPs_Res_P005UnRe.csv", sep=",", row.names=FALSE)
# 数百MB, 開けるのに時間かかる.
nrow(resOrdered_P005BaUn)
nrow(resOrdered_P005UnRe)


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
myAnnotation_BaUn <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                             analysis.type = "differential", design = design, 
                             contrasts = TRUE, cont.matrix = contMatrix, 
                             coef = "Ba - Un",
                             annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19"))
str(myAnnotation_BaUn)

myAnnotation_UnRe <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                                  analysis.type = "differential", design = design, 
                                  contrasts = TRUE, cont.matrix = contMatrix, 
                                  coef = "Un - Re",
                                  annotation=c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b2.hg19"))
str(myAnnotation_BaUn)

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
         what = "Beta", arraytype = "EPIC", genome = "hg19")

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

#### Gene Ontology Testing ####
# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
sigCpGs <- DMPs$Name[DMPs$P.Value<0.05]
# First 10 & Total number of significant CpGs at above threshold
sigCpGs[1:10]
length(sigCpGs)

# Get all the CpG sites used in the analysis to form the background
# & Total number of CpG sites tested
all <- DMPs$Name
length(all)

par(mfrow=c(1,1))
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all, plot.bias=TRUE)
# (Note that the warning regarding multiple symbols will always
# be displayed as there are genes that have more than one alias,
# however it is not a cause for concern.)

#. Top 10 GO Categories ----
t10GO <- topGSA(gst, number=10)
# write.csv(t10GO, "122120_go.csv")
