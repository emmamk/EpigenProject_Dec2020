#### Calculate M & beta-values ####
source("000_util.R")
# source("001_data_download.R")
# source("002_sort_pdata.R")
# source("003_quality_control.R")


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
    print(head(mVals[[i]][1:3,1:5]))
    cat("\n")
    
    # beta-values(getbeta(minfi)) ----
    # are easy to interpret: better for displaying data
    bVals[[i]] <- getBeta(mSetSqFlt[[i]])
    cat("beta-values:", "\n")
    print(head(bVals[[i]][1:3,1:5]))
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
    
    timepoint <- factor(pD[[i]]$timepoint)  # factor of interest
    individual <- factor(pD[[i]]$sample)  # individual effect that we need to account for
    
    # design matrix
    design <- model.matrix(~0+timepoint+individual, data = pD[[i]])
    cat("design[1:5,]:", "\n")
    print(design[1:5,])
    colnames(design) <- c(levels(timepoint), levels(individual)[-1])
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
    cat("summary(decideTests(fit2):", "\n")
    print(summary(decideTests(fit2[[i]])))
    cat("\n")
    # before - after
    # Down                0
    # NotSig         371592
    # Up                  0
    
    cat("summary(decideTests(fit2) adjP < 0.2:", "\n")
    print(summary(decideTests(fit2[[i]], p.value = 0.2)))
    cat("\n")
    # before - after
    # Down            12869
    # NotSig         348876
    # Up               9847
    }



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
    annotations[[i]] <- annotation[match(rownames(mVals[[i]]),annotation$Name), c(1:4,12:19,24:ncol(annotation))]
    names(annotations) <- c("ann450kSub", "annEPICSub")
    }



#### Get the Table of Results ####
# for the first contrast: coef  =  1
# (DMPs  =  differentially methylated probes )
# (To order by p-value, the user can specify topTable(sort.by = "p").
# DefultはB-statistic("B"). 但しmost casesでp値とB-statisticはidentical.)
# num = Inf はデータを全て摘出. Cut Offを設けていない.
DMPs <- list()
target_results <- list()
for (i in seq_along(fit2)) {
    annotation <- annotations[[i]]
    dmps <- topTable(fit2[[i]], num = Inf, coef = 1, genelist = annotation, sort.by = "p")
    DMPs[[i]] <- dmps
    
    cat(gseid[[i]], "\n")
    print(head(dmps, 3))
    cat("\n")
    #                   Regulatory_Feature_Group  DHS      logFC   AveExpr         t      P.Value adj.P.Val        B
    # cg11435841 Unclassified_Cell_type_specific TRUE  0.4434252 -1.119428  7.133892 7.253454e-07 0.1064706 5.480562
    # cg24968721                                       0.3575819  1.652057  7.076052 8.154710e-07 0.1064706 5.386572
    # cg14414873                                      -0.3506476  1.232118 -6.743198 1.614033e-06 0.1064706 4.834385


    ##### Explore Results #####
    cat("summary(dmps):", "\n")
    print(summary(dmps))
    cat("dmps$adj.P.Val < 0.2:", sum(dmps$adj.P.Val < 0.2, na.rm = TRUE), "\n")
    
    target_results[[i]] <- subset(dmps, adj.P.Val < 0.2)
    target_res <- target_results[[i]]
    print(head(target_res[order(target_res$logFC),], 3))
    print(head(target_res[order(-target_res$logFC),], 3))
    
    res_ordered <- target_res[order(target_res$adj.P.Val),]
    head(select(res_ordered, adj.P.Val))
    #            adj.P.Val
    # cg11435841 0.1064706
    # cg24968721 0.1064706
    # cg14414873 0.1064706
    
    # write.table(dmps, file = paste0("dmps_all_082023", gseid[i], ".csv"), sep = ",", row.names = FALSE)
    # write.table(res_ordered, file = paste0("dmps_adjp02_ordered_082023", gseid[i], ".csv"), sep = ",", row.names = FALSE)
    # 数百MB, 開けるのに時間かかる.
    
    
    #. Plot Top 4 Most Sig DM CpGs ----
    # 出力したいデータ  =  第一引数
    par(mfrow = c(2,2))
    sapply(rownames(res_ordered)[1:4],
           function(cpg){
               plotCpg(bVals, cpg = cpg, pheno = pD[[i]]$timepoint, ylab = "Beta values")
           }
           )
    }


#### Additional analyses ####
#. Gene Ontology Testing ####
# Get the significant CpG sites at less than FDR < 0.2
# 08/20/21: GSE114763とのOverlapsを抽出
head(DMPs)
GSE60655Sig <- subset(DMPs, adj.P.Val < 0.2)
GSE60655_P <- subset(GSE60655Sig, logFC > 0)  # nrow(GSE60655_P)  # 9847
GSE60655_N <- subset(GSE60655Sig, logFC < 0)  # nrow(GSE60655_N)  # 12869

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

par(mfrow = c(1,2))
gst_P <- gometh(sig.cpg = OLCpGs_P, all.cpg = all, plot.bias = TRUE)
gst_N <- gometh(sig.cpg = OLCpGs_N, all.cpg = all, plot.bias = TRUE)
# warning regarding multiple symbols will always be displayed as
# there are genes that have more than one alias; not a cause for concern.

#. . Top 10 GO Categories ----
t10GO_P <- topGSA(gst_P, num = Inf)
t10GO_P <- t10GO_P[t10GO_P$P.DE < 0.05,]
t10GO_N <- topGSA(gst_N, num = Inf)
t10GO_N <- t10GO_N[t10GO_N$P.DE < 0.05,]
length(t10GO_P$TERM)
length(t10GO_N$TERM)

t10GO_P$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO_P$TERM)]
t10GO_N$TERM[grep("blood pressure|cardio vascular", t10GO_N$TERM)]

t10GO <- rbind(t10GO_P, t10GO_N)
write.csv(t10GO, "go.csv")
