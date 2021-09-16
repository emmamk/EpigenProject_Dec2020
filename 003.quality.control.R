#### Quality Control ####
# ref1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3669124/
# ref2.: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4622000/

# setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/')
# source("000.util.R")
# source("001.data.prep.R")




for (i in seq_along(gseid)) {
    cat("#####", gsedir[[i]], "####", "\n")
    target.mlset <- mlset[[i]]
    
    #### removing 65 SNPs ####ß
    cat("removing 65 SNPs...", "\n")
    cat("dim(mlset) before removing SNPs:", dim(target.mlset), "\n")  # 485512  14
    
    snps <- featureNames(target.mlset)[grep("rs[0-9]+", featureNames(target.mlset))]
    keep <- !(featureNames(target.mlset) %in% snps)
    
    print(table(keep))
    # keep
    # FALSE   TRUE 
    # 65 485512
    
    target.mlset <- target.mlset[keep,]
    cat("dim(mlset) after removing SNPs:", dim(target.mlset), "\n\n")  # 485512  14
    
    
    #### removing probes on sex chromosomes (skipped: subjects are all male) ####
    
    
    #### removing probes by detection p-val #### 
    cat("checking detection p values...", "\n")
    detP <- target.mlset@assayData$pvals
    # length(detP[detP > 0.01])  # 2062
    # max(colMeans(detP))  # 0.0001693573
    
    #. removing detection p value > 0.01 in > 5% of the samples ----
    cat("table(rowSums(detP < 0.01) > 0.95*ncol(mlset))", "\n")
    print(table(rowSums(detP < 0.01) > 0.95*ncol(target.mlset)))  # 'FALSE' to be removed
    # FALSE   TRUE 
    #   978 484534
    
    keep <- rowSums(detP < 0.01) > 0.95*ncol(target.mlset)
    mlset[[i]] <- target.mlset[keep,]
    cat("dim(mlset) after removing poor quality probes:", dim(mlset[[i]]), "\n\n")  # 484534  14
    
}



#### removing outlier (b & M vals) ####
for (i in seq_along(gseid)) {
    cat("#####", gsedir[[i]], "####", "\n")
    cat("detecting outliers...", "\n")
    mVals.pre <- estimateM(as(mlset[[i]], "MethyLumiM"), returnType = "matrix")
    bVals.pre <- mlset[[i]]@assayData$betas
    
    par(mfrow=c(1,2))
    outlier.mvals <- detectOutlier(mVals.pre, ifPlot = TRUE)  # none
    outlier.bvals <- detectOutlier(bVals.pre, ifPlot = TRUE)  # none
    cat("check the destances between samples in the plot.", "\n\n")
    # S1.baseline in GSE114763 is detected as an outlier
    
    if (i != length(gsedir)) { next }
    
    cat("removing S1.baseline from GSE114763...", "\n")
    mlset[[gse]] <- mlset[[gse]][,colnames(mlset[[gse]]) != "S1.baseline"]
    cat("dim(mlset) after removing outlier:", dim(mlset[[gse]]))  # 451492       15
    pD[[gse]] <- pD[[gse]][rownames(pD[[gse]]) != "S1.baseline",]
    
    par(mfrow=c(1,1))
}


cat("\n\n")