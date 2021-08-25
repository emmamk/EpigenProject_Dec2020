#### Quality Control ####
source("000_util.R")
# source("001_data_download.R")
# source("002_sort_pdata.R")


#### Detecting p-values ####
detP <- list()
pal <- brewer.pal(8,"Dark2")
if (!dir.exists("plot")) {
    dir.create("plot")
}

png("plot/detectionp.png", width = 1280, height = 800)
par(mfrow = c(1,2))
for (i in seq_along(rgSet)) {
    cat(gseid[i], "\n")
    pdata <- pData(rgSet[[i]])
    detP[[i]] <- detectionP(rgSet[[i]])
    detection_p <- detP[[i]]
    
    cat("max(colMeans(detection_p)):", max(colMeans(detection_p)), "\n")  # 0.002766066
    cat("detP[1:3,1:5]:", "\n")
    print(detection_p[1:3,1:5])
    cat("\n\n")
    
    #. Plot detP (barplot) ----
    # examine mean detection p-values across all samples to identify any failed samples
    barplot(colMeans(detection_p), col = pal[factor(pdata$timepoint)], las = 2,
        cex.axis = 0.8, main = paste0(gseid[i], " detP"), ylab = "Mean detection p-values")
    abline(h = 0.05, col = "red")
    legend("topleft", legend = levels(factor(pdata$timepoint)), fill = pal, bg = "white", bty = "n") # cex = 0.5
    }
dev.off()

# memo ----
# GSE60655 
# max(colMeans(detection_p)): 0.0003285641 
# detP[1:3,1:5]: 
# S1.before      S1.after    S4.before      S4.after    S5.before
# cg00050873  0.000000e+00  0.000000e+00  0.00000e+00  0.000000e+00 0.000000e+00
# cg00212031  0.000000e+00  0.000000e+00  0.00000e+00  0.000000e+00 0.000000e+00
# cg00213748 3.333207e-103 2.922029e-120 3.99785e-104 3.585549e-146 2.210492e-85
# 
# 
# GSE114763 
# max(colMeans(detection_p)): 0.0007908372 
# detP[1:3,1:5]: 
# S1.before     S1.after S2.before S2.after     S3.before
# cg18478105         0  0.00000e+00         0        0  0.000000e+00
# cg09835024         0 7.22592e-279         0        0 6.368765e-263
# cg14361672         0  0.00000e+00         0        0  0.000000e+00



#### qcReport(minfi) ####
# Ex. of another useful quality control plot
# 450K用? 850Kは不明
# library(minfi)
# for (i in seq_along(rgSet)) {
#     name <- gseid[i]
#     pdata <- pData(rgSet[[i]])
#     qcReport(rgSet[[i]], sampNames = pdata$ID, sampGroups = pdata$timepoint, 
#              pdf = paste0("plot/qcReport_", name, ".pdf"))
#     dev.off()
#     }

#### Remove Poor Quality Samples ####
for (i in seq_along(detP)) {
    cat(gseid[i], "\n")
    cat("max(colMeans(detP):", max(colMeans(detP[[i]])), "\n")
    keep <- colMeans(detP[[i]]) < 0.05
    cat("samples to be kept (table(keep)):", table(keep), "\n")
    # TRUE 
    # 34
    
    # removing low detP data from rgSet & pData ----
    rgSet[[i]] <- rgSet[[i]][,keep]
    cat("dim(rgSet):", dim(rgSet[[i]]), "\n")
    pD[[i]] <- pD[[i]][keep,]
    cat("dim(pD):", dim(pD[[i]]), "\n")
    
    # removing above from detection p-value table ----
    detP[[i]] <- detP[[i]][,keep]
    cat("dim(detP):", dim(detP[[i]]), "\n\n")  # 485512     34
}

# memo: 
# GSE60655 
# max(colMeans(detP): 0.0003285641 
# samples to be kept (table(keep)): 14 
# dim(rgSet): 622399 14 
# dim(pD): 14 6 
# dim(detP): 485512 14 
#     
# GSE114763 
# max(colMeans(detP): 0.0007908372 
# samples to be kept (table(keep)): 16 
# dim(rgSet): 1051943 16 
# dim(pD): 16 6 
# dim(detP): 866238 16



#### Normalisation ####
# normalize the data; this results in a GenomicRatioSet object
# choose optimal normalization method suitable for the data
# ("Warning messages:" regarding "only one sex is present..." can be ignored.)
mSetSq <- list()
mSetRaw <- list()
for (i in seq_along(rgSet)) {
    cat(gseid[i], "\n")
    #. normalizing ----
    mSetSq[[i]] <- preprocessQuantile(rgSet[[i]])
    # create a MethylSet object from the raw data for later plotting
    mSetRaw[[i]] <- preprocessRaw(rgSet[[i]])
    
    #. densityplot: Before & After Normalisation ----
    png(paste0("plot/densityplot_", gseid[i], ".png"), width = 1280, height = 800)
    par(mfrow = c(1,2))
    pdata <- pD[[i]]
    # Before
    densityPlot(rgSet[[i]], sampGroups = pdata$timepoint, main = paste0(gseid[i], " Raw"), legend = FALSE)
    legend("top", legend  =  levels(factor(pdata$timepoint)),
           text.col = brewer.pal(8,"Dark2"), cex = 0.8, bty = "n")
    
    # After
    densityPlot(getBeta(mSetSq[[i]]), sampGroups = pdata$timepoint, main = "Normalized", legend = FALSE)
    legend("top", legend  =  levels(factor(pdata$timepoint)),
           text.col = brewer.pal(8,"Dark2"), cex = 0.8, bty = "n")
    dev.off()
    cat("\n\n")
}



#### Data exploration ####
#### MDS plot ####
# Multi-dimensional scaling (MDS) plots: to look at largest sources of variation
for (i in seq_along(mSetSq)) {
    cat(gseid[i], "\n")
    png(paste0("plot/mds_", gseid[i], ".png"))
    plotMDS(getM(mSetSq[[i]]), top = 1000, gene.selection = "common", 
            col = pal[factor(pD[[i]]$timepoint)])
    legend("topleft", legend = levels(factor(pD[[i]]$timepoint)), text.col = pal,
           fill = pal, bg = "white", cex = 0.8, bty = "n")
    dev.off()
    
    # plotMDS(getM(mSetSq), top = 1000, gene.selection = "common",  
    #         col = pal[factor(pD$sample)], cex = 0.6)
    # legend("top", legend = levels(factor(pD$sample)), text.col = pal,
    #      fill = pal, bg = "white", cex = 0.5, bty = "n")
    
    #. . Examine Higher dimensions ----
    png(paste0("plot/mds2_", gseid[i], ".png"), width = 800)
    par(mfrow = c(1,3))
    plotMDS(getM(mSetSq[[i]]), top = 1000, gene.selection = "common", 
            col = pal[factor(pD[[i]]$timepoint)], dim = c(1,3), cex = 1.2)
    legend("topleft", legend = levels(factor(pD[[i]]$timepoint)), text.col = pal,
           fill = pal, bg = "white", bty = "n")
    plotMDS(getM(mSetSq[[i]]), top = 1000, gene.selection = "common", 
            col = pal[factor(pD[[i]]$timepoint)], dim = c(2,3), cex = 1.2)
    plotMDS(getM(mSetSq[[i]]), top = 1000, gene.selection = "common", 
            col = pal[factor(pD[[i]]$timepoint)], dim = c(3,4), cex = 1.2)
    dev.off()
    }



#### Filtering ####
#. Poor Performing Probes ----
# Poor performing probes(unreliable signal) are generally filtered out
# prior to differential methylation analysis.
# i.e. fewer statistical tests
mSetSqFlt <- list()
for (i in seq_along(detP)) {
    cat(gseid[i], "\n")
    # ensure probes are in the same order in the mSetSq and detP objects
    detP[[i]] <- detP[[i]][match(featureNames(mSetSq[[i]]), rownames(detP[[i]])),]
    
    # remove any probes that have failed in one or more samples
    # ex.detPのcolに1つNAがあるとrowSums(detP)は9、一方ncol(mSetSq)は常に10.
    # そのように左右の数が合わない物を除外する.
    keep <- rowSums(detP[[i]] < 0.01) == ncol(mSetSq[[i]])
    cat("probes to be kept (table(keep)):", table(keep), "\n")
    mSetSqFlt[[i]] <- mSetSq[[i]][keep,]
    cat("after removing uncomplete samples (dim(mSetSqFlt)):", dim(mSetSqFlt[[i]]), "\n")
    
    #. Probes with SNPs at CpG Site ----
    # (common SNPs may affect the CpG)
    mSetSqFlt[[i]] <- dropLociWithSnps(mSetSqFlt[[i]])
    cat("after removing common SNPs (dim(mSetSqFlt)):", dim(mSetSqFlt[[i]]), "\n\n")
    }

# GSE60655 
# probes to be kept (table(keep)): 1290 484222 
# after removing uncomplete samples (dim(mSetSqFlt)): 484222 14 
# after removing common SNPs (dim(mSetSqFlt)): 467244 14 
# 
# GSE114763 
# probes to be kept (table(keep)): 6504 859355 
# after removing uncomplete samples (dim(mSetSqFlt)): 859355 16 
# after removing common SNPs (dim(mSetSqFlt)): 829814 16 



#### Cross Reactive Probes ####
source("004_removing_cross_reactive_probes_450k.R")
source("004_removing_cross_reactive_probes_EPIC.R")



#### Re-examine Data with MDS Plots ####
# Once the data has been filtered and normalised, it is often useful to re-examine
# the MDS plots to see if the relationship between the samples has changed
for (i in seq_along(mSetSqFlt)) {
    cat(gseid[i], "\n")
    
    png(paste0("plot/mds_afterqc_", gseid[i], ".png"))
    par(mfrow = c(1,1))
    plotMDS(getM(mSetSqFlt[[i]]), top = 1000, gene.selection = "common", 
            col = pal[factor(pD[[i]]$timepoint)], main = paste0(gseid[i], " After Filtering"))
    legend("right", legend = levels(factor(pD[[i]]$timepoint)), text.col = pal,
           cex = 0.8, fill = pal, bg = "white",  bty = "n")
    dev.off()
    
    png(paste0("plot/mds_bf.af_qcfiltering", gseid[i], ".png"), width = 800)
    par(mfrow = c(1,2))
    # Before
    plotMDS(getM(mSetSq[[i]]), top = 1000, gene.selection = "common", main = paste0(gseid[i], " Before"), 
            col = pal[factor(pD[[i]]$timepoint)])
    legend("topleft", legend = levels(factor(pD[[i]]$timepoint)), text.col = pal,
           cex = 0.8, bg = "white", fill = pal, bty = "n")
    # After
    plotMDS(getM(mSetSqFlt[[i]]), top = 1000, gene.selection = "common", main = "After",
            col = pal[factor(pD[[i]]$timepoint)])
    dev.off()
}
