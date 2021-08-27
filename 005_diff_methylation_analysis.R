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
    cat("####", gsedir[i], "####", "\n")
    gse <- gseid[i]
    # M-values(getM(minfi)) ----
    mVals[[gse]] <- getM(mSetSqFlt[[i]])
    cat("m-values:", "\n")
    print(mVals[[i]][1:3,1:5])
    cat("\n")
    
    # beta-values(getbeta(minfi)) ----
    # are easy to interpret: better for displaying data
    bVals[[gse]] <- getBeta(mSetSqFlt[[i]])
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
#### GSE60655_End #### 
# m-values: 
#             S1.before  S1.after  S4.before   S4.after S5.before
# cg13869341  2.5053762  2.566522  2.5152105  2.7193212  2.359755
# cg14008030  0.9984332  1.082772  0.8438259  0.7416271  1.054326
# cg12045430 -1.7682828 -1.151728 -2.3306506 -2.4950569 -2.148532
# 
# beta-values: 
#            S1.before  S1.after S4.before  S4.after S5.before
# cg13869341 0.8502540 0.8555705 0.8511198 0.8681718 0.8369432
# cg14008030 0.6664253 0.6792921 0.6421931 0.6257570 0.6749814
# cg12045430 0.2269384 0.3103857 0.1658286 0.1506590 0.1840345
# 
# 
#### GSE114763_Res #### 
# m-values: 
#            S1.after S2.before  S2.after S3.before  S3.after
# cg14817997 1.400965 0.8085359 0.3092162 0.5276493 0.7237121
# cg26928153 2.460918 2.7322902 2.9874624 2.4053718 1.8832018
# cg16269199 1.083570 1.0241727 1.2090697 0.9312085 0.7509296
# 
# beta-values: 
#             S1.after S2.before  S2.after S3.before  S3.after
# cg14817997 0.7253337 0.6365530 0.5533789 0.5904289 0.6228444
# cg26928153 0.8462879 0.8691972 0.8880277 0.8412124 0.7867317
# cg16269199 0.6794126 0.6703796 0.6980572 0.6559877 0.6272658



#### Probe-wise Differential Methylation Analysis ####
#. Create a Design Matrix ----
# library(DMRcate)
# library(limma)
fit2 <- list()
for (i in seq_along(mVals)) {
    cat("####", gsedir[i], "####", "\n")
    gse <- gseid[i]
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
    fit2[[gse]] <- contrasts.fit(fit, contMatrix)
    fit2[[gse]] <- eBayes(fit2[[i]])
    
    #. Numbers of DM CpGs at FDR < 0.05 ----
    cat("####", gsedir[i], "####", "\n")
    cat("FDR < 0.2:", "\n")
    print(summary(decideTests(fit2[[i]], p.value = 0.2)))
    cat("\n")
    
    cat("FDR < 0.05", "\n")
    print(summary(decideTests(fit2[[i]])))  # default
    cat("\n")
    
    cat("p < 0.05:", "\n")
    print(summary(decideTests(fit2[[i]], adjust.method = "none", p.value = 0.05)))
    cat("\n")
    
    cat("p < 0.01:", "\n")
    print(summary(decideTests(fit2[[i]], adjust.method = "none", p.value = 0.01)))
    cat("\n\n")
    }

# memo:
#### GSE60655_End ####
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
# before - after
# Down            16459
# NotSig         352598
# Up              16163
# 
# p < 0.01: 
# before - after
# Down             2911
# NotSig         379033
# Up               3276

#### GSE114763_Res (Base vs Loading) ####
# FDR < 0.2: 
# before - after
# Down                0
# NotSig         777574
# Up                  0
# 
# FDR < 0.05 
# before - after
# Down                0
# NotSig         777574
# Up                  0
# 
# p < 0.05: 
# before - after
# Down            19709
# NotSig         737501
# Up              20364
# 
# p < 0.01: 
# before - after
# Down             3297
# NotSig         770604
# Up               3673



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
    gse <- gseid[i]
    annotation <- annotations[[i]]
    dmps <- topTable(fit2[[i]], num = Inf, coef = 1, genelist = annotation, sort.by = "p")
    DMPs[[gse]] <- dmps
    
    cat("####", gsedir[i], "####", "\n")
    print(dmps[1:3, c("logFC", "P.Value", "adj.P.Val")])
    cat("\n")
    }

#### GSE60655_End #### 
#                 logFC      P.Value adj.P.Val
# cg23093333  0.6945042 3.182226e-06 0.5720352
# cg02052895 -0.7176025 8.672024e-06 0.5720352
# cg04778337  0.8897548 1.137924e-05 0.5720352

#### GSE114763_Res #### 
#                 logFC      P.Value adj.P.Val
# cg10474429  0.6777285 2.962484e-06 0.9256469
# cg05588658  0.8273751 3.695408e-06 0.9256469
# cg14838474 -0.5963328 7.618736e-06 0.9256469



##### Explore Results #####
resultsSig <- list()
for (i in seq_along(DMPs)) {
    dmps <- DMPs[[i]]
    cat("#### results summary:", gsedir[i], "####", "\n")
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

    # write.table(dmps, file = paste0("res.all.0825.", gseid[i], ".csv"), sep = ",", row.names = FALSE)
    # write.table(res_p05, file = paste0("res.p05.0825.", gseid[i], ".csv"), sep = ",", row.names = FALSE)
    # write.table(res_p01, file = paste0("res.p01.0825.", gseid[i], ".csv"), sep = ",", row.names = FALSE)
    }

# memo:
# #### results summary: GSE60655_End #### 
# adj.p < 0.2: 0 pos: 0 neg: 0 
# p < 0.05: 40043 pos: 19875 neg: 20168 
# p < 0.01: 7674 pos: 4081 neg: 3593 
# 
#### results summary: GSE114763_Res ####
# adj.p < 0.2: 0 pos: 0 neg: 0 
# p < 0.05: 42951 pos: 21640 neg: 21311 
# p < 0.01: 7476 pos: 3915 neg: 3561 


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

