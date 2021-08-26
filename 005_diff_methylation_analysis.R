#### Scripts to get results table ####
source("000_util.R")
source("001_data_download.R")
source("002_sort_pdata.R")
source("003_quality_control.R")


#### Get M/beta-values ####
# nicer statistical properties: better for use in statistical analysis
# library(minfi)
mVals <- list()
bVals <- list()
for (i in seq_along(mSetSqFlt)) {
    cat(gseid[i], "\n")
    # M-values(getM(minfi)) ----
    mVals[[i]] <- getM(mSetSqFlt[[i]])
    cat("m-values:", "\n")
    print(mVals[[i]][1:3,1:5])
    cat("\n")
    
    # beta-values(getbeta(minfi)) ----
    # are easy to interpret: better for displaying data
    bVals[[i]] <- getBeta(mSetSqFlt[[i]])
    cat("beta-values:", "\n")
    print(bVals[[i]][1:3,1:5])
    cat("\n\n")
    
    #. Visualise M & beta-values ----
    png(paste0("plot/m.beta.values_", gseid[i], ".png"))
    par(mfrow = c(1,2))
    densityPlot(mVals[[i]], sampGroups = pD[[i]]$timepoint, main = paste0(gseid[i], " M-values"), 
                legend = FALSE, xlab = "M values")
    densityPlot(bVals[[i]], sampGroups = pD[[i]]$timepoint, main = "Beta values", 
                legend = FALSE, xlab = "Beta values")
    legend("top", legend  =  levels(factor(pD[[i]]$timepoint)),
           text.col = brewer.pal(8,"Dark2"), bty = "n")
    dev.off()
}

# memo:
# GSE60655 
# m-values: 
#            S1.before S1.after S4.before S4.after S5.before
# cg13869341  2.505376 2.566522  2.515211 2.719321  2.359755
# cg24669183  2.182828 2.000885  1.839788 2.179010  1.814513
# cg15560884  1.384137 1.322324  1.414610 1.109830  1.292936
# 
# beta-values: 
#            S1.before  S1.after S4.before  S4.after S5.before
# cg13869341 0.8502540 0.8555705 0.8511198 0.8681718 0.8369432
# cg24669183 0.8195085 0.8000982 0.7816391 0.8191167 0.7786342
# cg15560884 0.7230038 0.7143417 0.7272139 0.6833641 0.7101670
# 
# 
# GSE114763 
# m-values: 
#            S1.before S1.after S2.before S2.after S3.before
# cg26928153 2.1469000 2.440795  2.711936 2.966014 2.3866463
# cg16269199 0.7596354 1.063781  1.005536 1.192596 0.9129245
# cg13869341 1.4357869 1.607879  1.633010 1.736525 1.5704599
# 
# beta-values: 
#            S1.before  S1.after S2.before  S2.after S3.before
# cg26928153 0.8157956 0.8444647 0.8675848 0.8865408 0.8394710
# cg16269199 0.6286756 0.6764176 0.6675188 0.6956449 0.6531221
# cg13869341 0.7301161 0.7529665 0.7561924 0.7691766 0.7481104



#### Probe-wise Differential Methylation Analysis ####
#. Create a Design Matrix ----
# library(DMRcate)
# library(limma)
fit2 <- list()
for (i in seq_along(mVals)) {
    cat(gsedir[i], "\n")
    timepoint <- factor(pD[[i]]$timepoint)  # factor of interest
    individual <- factor(pD[[i]]$sample)  # individual effect that we need to account for
    
    # design matrix
    design <- model.matrix(~0+timepoint+individual, data = pD[[i]])
    cat("design[1:5,]:", "\n")
    print(design[1:5,])
    cat("\n")
    
    colnames(design) <- c(levels(timepoint), levels(individual)[-1])
    cat("renamed design[1:5,]:", "\n")
    print(design[1:5,])
    cat("\n")
    
    #. Fit the Linear Model ----
    # lmFit: Linear Model For Series Of Arrays
    fit <- lmFit(mVals[[i]], design)
    cat("mVals[1:3, 1:5]:", "\n")
    print(mVals[[i]][1:3, 1:5])
    cat("\n")
    #                 S1.b      S1.a      S4.b      S4.a      S5.b
    # cg13869341  2.578282  2.641601  2.588604  2.799171  2.427225
    # cg24669183  2.246692  2.056624  1.883342  2.237400  1.863289
    # cg15560884  1.412289  1.346616  1.446557  1.122430  1.313780
    
    #. Contrast Matrix ----
    # for specific comparisons
    contMatrix <- makeContrasts(before-after, levels = design)
    
    #. Fit the Contrasts ----
    fit2[[i]] <- contrasts.fit(fit, contMatrix)
    fit2[[i]] <- eBayes(fit2[[i]])
    
    #. Numbers of DM CpGs at FDR < 0.05 ----
    cat(gsedir[i], "\n")
    cat("FDR < 0.2:", "\n")
    print(summary(decideTests(fit2[[i]], p.value = 0.2)))
    cat("\n")
    # before - after
    # Down            12869
    # NotSig         348876
    # Up               9847
    
    cat("FDR < 0.05", "\n")
    print(summary(decideTests(fit2[[i]])))  # default
    cat("\n")
    # before - after
    # Down                0
    # NotSig         371592
    # Up                  0
    
    
    cat("p < 0.05:", "\n")
    print(summary(decideTests(fit2[[i]], adjust.method = "none", p.value = 0.05)))
    cat("\n")
    # before - after
    # Down            12869
    # NotSig         348876
    # Up               9847
    
    cat("p < 0.01:", "\n")
    print(summary(decideTests(fit2[[i]], adjust.method = "none", p.value = 0.01)))
    cat("\n\n")
    # before - after
    # Down            12869
    # NotSig         348876
    # Up               9847
    }


# memo:
#### GSE60655_End 
# FDR < 0.2: 
#     before - after
# Down                0
# NotSig         385220
# Up                  0
# 
# FDR < 0.05 
# before - after
# Down                0
# NotSig         385220
# Up                  0
# 
# p < 0.05: 
#     before - after
# Down            16459
# NotSig         352598
# Up              16163
# 
# p < 0.01: 
#     before - after
# Down             2911
# NotSig         379033
# Up               3276

#### GSE114763_Res 
# FDR < 0.2: 
#     before - after
# Down                0
# NotSig         776698
# Up                  0
# 
# FDR < 0.05 
# before - after
# Down                0
# NotSig         776698
# Up                  0
# 
# p < 0.05: 
#     before - after
# Down            18058
# NotSig         741895
# Up              16745
# 
# p < 0.01: 
#     before - after
# Down             2766
# NotSig         771419
# Up               2513



#### Get Annotation Data ####
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(IlluminaHumanMethylation450kmanifest)
# library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
# library(IlluminaHumanMethylationEPICmanifest)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
annotations <- list(ann450k, annEPIC)

for (i in seq_along(annotations)) {
    annotation <- annotations[[i]]
    annotations[[i]] <- 
        annotation[match(rownames(mVals[[i]]),annotation$Name), c(1:4,12:19,24:ncol(annotation))]
    names(annotations) <- c("ann450kSub", "annEPICSub")
    }



#### Get the Table of Results ####
# for the first contrast: coef  =  1
# (DMPs  =  differentially methylated probes )
# (To order by p-value, the user can specify topTable(sort.by = "p").
# DefultはB-statistic("B"). 但しmost casesでp値とB-statisticはidentical.)
# num = Inf はデータを全て摘出. Cut Offを設けていない.
DMPs <- list()
for (i in seq_along(fit2)) {
    annotation <- annotations[[i]]
    dmps <- topTable(fit2[[i]], num = Inf, coef = 1, genelist = annotation, sort.by = "p")
    DMPs[[i]] <- dmps
    
    cat(gsedir[i], "\n")
    print(dmps[1:3, c("logFC", "P.Value", "adj.P.Val")])
    cat("\n")
    #                   Regulatory_Feature_Group  DHS      logFC   AveExpr         t      P.Value adj.P.Val        B
    # cg11435841 Unclassified_Cell_type_specific TRUE  0.4434252 -1.119428  7.133892 7.253454e-07 0.1064706 5.480562
    # cg24968721                                       0.3575819  1.652057  7.076052 8.154710e-07 0.1064706 5.386572
    # cg14414873                                      -0.3506476  1.232118 -6.743198 1.614033e-06 0.1064706 4.834385
    }

##### Explore Results #####
resultsSig <- list()
for (i in seq_along(DMPs)) {
    dmps <- DMPs[[i]]
    cat("results summary:", gsedir[i], "\n")
    cat("adj.p < 0.2:", sum(dmps$adj.P.Val < 0.2, na.rm = TRUE),
        "pos:", sum(dmps$adj.P.Val < 0.2 & dmps$logFC > 0, na.rm = TRUE),
        "neg:", sum(dmps$adj.P.Val < 0.2 & dmps$logFC < 0, na.rm = TRUE), "\n")
    cat("p < 0.05:", sum(dmps$P.Value < 0.05, na.rm = TRUE),
        "pos:", sum(dmps$P.Value < 0.05 & dmps$logFC > 0, na.rm = TRUE),
        "neg:", sum(dmps$P.Value < 0.05 & dmps$logFC < 0, na.rm = TRUE), "\n")
    cat("p < 0.01:", sum(dmps$P.Value < 0.01, na.rm = TRUE),
        "pos:", sum(dmps$P.Value < 0.01 & dmps$logFC > 0, na.rm = TRUE),
        "neg:", sum(dmps$P.Value < 0.01 & dmps$logFC < 0, na.rm = TRUE), "\n\n")
    
    resultsSig[[i]] <- subset(dmps, adj.P.Val < 0.2)
    resultsSig[[i+2]] <- res_p05 <- subset(dmps, P.Value < 0.05)
    resultsSig[[i+4]] <- res_p01 <- subset(dmps, P.Value < 0.01)
    names(resultsSig)[c(i, i+2, i+4)] <- paste0(rep(paste0("res", gseid[i]), 3), c("adj02", "p005", "p001"))

    write.table(dmps, file = paste0("res.all.0825.", gseid[i], ".csv"), sep = ",", row.names = FALSE)
    write.table(res_p05, file = paste0("res.p05.0825.", gseid[i], ".csv"), sep = ",", row.names = FALSE)
    write.table(res_p01, file = paste0("res.p01.0825.", gseid[i], ".csv"), sep = ",", row.names = FALSE)
    }

# memo:
# results summary: GSE60655_End 
# adj.p < 0.2: 0 pos: 0 neg: 0 
# p < 0.05: 32622 pos: 16163 neg: 16459 
# p < 0.01: 6187 pos: 3276 neg: 2911 
# 
# results summary: GSE114763_Res 
# adj.p < 0.2: 0 pos: 0 neg: 0 
# p < 0.05: 34803 pos: 16745 neg: 18058 
# p < 0.01: 5279 pos: 2513 neg: 2766


# check if named correctly
names(resultsSig)
# [1] "resGSE60655adj02"  "resGSE114763adj02" "resGSE60655p005"   "resGSE114763p005"  "resGSE60655p001"  
# [6] "resGSE114763p001"




#### Plot Top 4 Most Sig DM CpGs (optional) ####
# cpg: spot to plot
# par(mfrow = c(2,2), oma = c(0,0,4,0))
# for (i in seq_along(gseid)) {    
#     res <- resultsSig[[i]]
#     res_ordered <- res[order(res$adj.P.Val),]  # just in case
#     sapply(rownames(res_ordered)[1:4],
#            function(cpg){
#                print(
#                    plotCpg(bVals[[i]], cpg = cpg, pheno = pD[[i]]$timepoint, ylab = "Beta values"),
#                    mtext(side = 3, line=1, outer=T, text = gsedir[i], cex=2)
#                )
#            }
#     )
# }
# par(mfrow = c(1,1), oma = c(0,0,4,0))

