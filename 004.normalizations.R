#### Normalizations ####
# project.dir <- "/Users/Emma/GitHub/EpigenProject_Dec2020/"
# setwd(project.dir)
# source("000.util.R")
# source("001.data.prep.R")
# source("002.pdata.sort.R")
# source("003.quality.control.R")


mlset.norm <- list()
combat <- list()
for (i in seq_along(gseid)) {
    cat("#####", gsedir[[i]], "####", "\n")
    pd <- pD[[i]]
    gse <- gseid[i]
    
    #### Color-bias adjustment ####
    cat("adjColorBias.quantile: adjusting color-bias...", "\n")
    mlset.clradj <- adjColorBias.quantile(mlset[[i]])
    par(mfrow = c(1,2))
    plotColorBias1D(mlset[[i]])
    plotColorBias1D(mlset.clradj)
    cat("\n")
    
    
    #### BMIQ (probe type bias adjustment) ####
    # normalization on β-values
    # library(wateRmelon)
    cat("BMIQ: normalization on β-values...", "\n")
    mlset.clradj.bmiq <- BMIQ(mlset.clradj)
    mlset.norm[[gse]] <- mlset.clradj.bmiq
    cat("comparing beta values density before & after BMIQ...", "\n\n")
    par(mfrow = c(1,2))
    plotDensity(mlset[[i]]@assayData$betas, legend = FALSE, main = "raw", col=pal[factor(pd$ID)])
    legend("topleft", legend=levels(factor(pd$ID)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")
    plotDensity(mlset.clradj.bmiq@assayData$betas, legend = FALSE, main = "BMIQed", col=pal[factor(pd$ID)])
    legend("topleft", legend=levels(factor(pd$ID)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")
    
    
    
    #### Data exploration ####
    #. MDS plot ----
    # Multi-dimensional scaling (MDS) plots: to look at largest sources of variation
    cat("exploring normalized data...", "\n")
    mVals <- estimateM(as(mlset.clradj.bmiq, "MethyLumiM"), returnType = "matrix")
    bVals <- mlset.clradj.bmiq@assayData$betas
    
    cat("MDS plots: checking PC1 & PC2...", "\n")
    par(mfrow=c(1,3))
    plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pd$timepoint)], cex=0.8)
    legend("bottomright", legend=levels(factor(pd$timepoint)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")
    
    plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pd$sample)], cex=0.8)
    legend("bottomright", legend=levels(factor(pd$sample)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")
    
    plotMDS(mVals, top=1000, gene.selection="common", col=pal[factor(pd$slide)], cex=0.8)
    legend("bottomright", legend=levels(factor(pd$slide)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")

        
    #. density plot: M & beta values ----
    cat("checking density of M & beta values...", "\n\n")
    par(mfrow=c(1,2))
    plotDensity(mVals, legend = FALSE, main = "M values", col=pal[factor(pd$timepoint)])
    legend("topright", legend=levels(factor(pd$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")
    plotDensity(bVals, legend = FALSE, main = "Beta values", col=pal[factor(pd$timepoint)])
    legend("topleft", legend=levels(factor(pd$timepoint)), text.col=pal, cex=0.8, bg="white", bty = "n")
    
    
    
    #### ComBat: removing batch effect ####
    cat("ComBat: removing batch effect (slide) ...", "\n")
    batch <- as.factor(pd$slide)
    modcombat <- model.matrix(~ 1, data = pd)
    combat.edata <- ComBat(dat = mVals, batch = batch, mod = modcombat, 
                           par.prior = TRUE, prior.plots = TRUE)
    combat[[gse]] <- combat.edata
    
    
    #. MDS: visualizing Combat()ed data ----
    par(mfrow=c(1,3))
    plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pd$timepoint)], cex=0.8)
    legend("bottomright", legend=levels(factor(pd$timepoint)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")
    
    plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pd$sample)], cex=0.8)
    legend("bottomright", legend=levels(factor(pd$sample)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")
    
    plotMDS(combat.edata, top=1000, gene.selection="common", col=pal[factor(pd$slide)], cex=0.8)
    legend("bottomright", legend=levels(factor(pd$slide)), text.col=pal, fill=pal, bg="white", cex=0.7, bty = "n")
    
}
    
    
cat("\n\n")