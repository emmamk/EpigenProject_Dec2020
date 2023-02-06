#### Probe-wise Differential Methylation Analysis ####
# project.dir <- "/Users/Emma/GitHub/EpigenProject_Dec2020/"
# setwd(project.dir)
# source("000.util.R")
# source("001.data.prep.R")
# source("002.pdata.sort.R")
# source("003.quality.control.R")
# source("004.normalizations.R")


cat("#### Differential Methylation Analysis Results Summary ####")
fit.list <- list()
for (i in seq_along(gseid)) {
    cat("#####", gsedir[[i]], "####", "\n")
    gse <- gseid[[i]]
    pd <- pD[[i]]

    #. Create a Design Matrix ----
    # ref: https://genomicsclass.github.io/book/pages/expressing_design_formula.html
    timepoint <- factor(pd$timepoint)  # factor of interest
    individual <- factor(pd$sample)  # adjustment
    design <- model.matrix(~ 0 + timepoint + individual, data = pd)
    colnames(design) <- c(levels(timepoint),levels(individual)[-1])
    
    #. Fit the Linear Model ----
    # lmFit: Linear Model For Series Of Arrays
    cat("fitting the data by lmFit...", "\n")
    fit <- lmFit(combat[[i]], design)
    
    cat("contrasting the model by before vs. after (eBayse)...", "\n")
    #. Contrast Matrix ----
    # for specific comparisons
    contMatrix <- makeContrasts(after-before, levels=design)  # ref: before
    # contMatrix
    
    #. Fit the Contrasts ----
    fit2 <- contrasts.fit(fit, contMatrix)
    fit2 <- eBayes(fit2)
    fit.list[[gse]] <- fit2
    
    #. Numbers of DM CpGs at FDR < 0.05 ----
    cat("summary(decideTests(fit2): Down, NotSig, Up", "\n")
    cat("FDR < 0.05:", summary(decideTests(fit2)), "\n")
    cat("FDR < 0.2:", summary(decideTests(fit2, p.value = 0.2)), "\n")
    cat("p < 0.05:", summary(decideTests(fit2, adjust.method = "none", p.value = 0.05)), "\n\n")
}    



#### Get the Table of Results ####
#. annotation ----
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

DMPs <- list()
DMPs005 <- list()
for (i in seq_along(gseid)) {
    cat("#####", gsedir[[i]], "####", "\n")
    gse <- gseid[[i]]
    ann450kSub <- ann450k[match(rownames(combat[[i]]), ann450k$Name), c(1:4,12:19,24:ncol(ann450k))]
    
    dmps <- topTable(fit.list[[i]], num=Inf, coef=1, genelist=ann450kSub, sort.by = "p")
    DMPs[[gse]] <- dmps
    
    print(dmps[1:3,c("logFC", "P.Value", "adj.P.Val")]); cat("\n")
    #                 logFC      P.Value adj.P.Val
    # cg08642731  0.6338006 6.213410e-06 0.4997675
    # cg02586025  0.5055961 2.447523e-05 0.4997675
    # cg12337011 -0.5370560 3.103759e-05 0.4997675
    
    dmps.sig <- dmps %>%
        select(logFC, P.Value, adj.P.Val) %>% 
        filter(P.Value < 0.05 & abs(logFC) > 1) %>% 
        nrow()  # 3
    cat("number of p < 0.05 & abs(logFC) > 1: ", dmps.sig , "\n\n")
    
    cat("summary(DMPs$logFC):", "\n")
    print(summary(dmps$logFC)); cat("\n")
    #      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
    # -1.222968 -0.097350  0.001802  0.004197  0.102325  1.169221
    
    
    # write.table(DMPs, file="DMPs.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv", sep=",", row.names=FALSE)
    dmps.005 <- subset(dmps, P.Value < 0.05)
    DMPs005[[gse]] <- dmps.005
    # write.table(dmps.005, file="DMPs.p005.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv", sep=",", row.names=FALSE)
}


#### Volcano plot ####
#### Volcano plot ####
# library(EnhancedVolcano)
# for (i in seq_along(DMPs)) {
#     dmps <- DMPs[[i]]
#     print(
#         EnhancedVolcano(dmps,
#                         lab = dmps$Name,
#                         x = 'logFC', FCcutoff = 0.5, #xlim = c(-2, 2),
#                         y = 'P.Value', pCutoff = 0.05, #ylim = c(0, 6),
#                         xlim = c(floor(min(dmps$logFC)), ceiling(max(abs(dmps$logFC)))),
#                         ylim = c(0, ceiling(max(-log10(dmps$P.Value)))),
#                         xlab = bquote(~Log[2]~ 'fold change'),
#                         ylab = bquote(~-Log[10]~italic(P)),  # ~adjusted
#                         pointSize = 2.0, 
#                         labSize = 3.0,
#                         col = c('black', 'red', 'orange', 'blue'),
#                         legendLabels = c('NS','Log2 FC','padj','padj & Log2 FC')) +
#             
#             labs(title = paste0(gsedir[[i]], " p < 0.05"), titleLabSize = 10,
#                  subtitle = "pre vs post", caption = "")
#     )
# }


#### Plot Top 4 Most Sig DM.CpGs ####
# cf. DMPs is already ordered by P.Value
# par(mfrow=c(2,2))
# sapply(rownames(dmps.005)[1:4], function(cpg){
#     plotCpg(bVals, cpg=cpg, pheno=pD$timepoint, ylab = "Beta values")
# })



cat("\n\n")