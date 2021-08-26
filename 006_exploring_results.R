#### Additional analyses ####
# source("000_util.R")
# source("001_data_download.R")
# source("002_sort_pdata.R")
# source("003_quality_control.R")
# source("005_diff_methylation_analysis.R")


#### Gene Ontology Testing ####
# Get the significant CpG sites
# FDR < 0.2 --- 0  (8/23/21)
# p < 0.01 ---   (8/25/21)

names(resultsSig)
# "resGSE60655adj02"  "resGSE114763adj02" "resGSE60655p005"   "resGSE114763p005" 
# "resGSE60655p001"  "resGSE114763p001"

GSE60655.sig <- resultsSig[["resGSE60655p005"]]
GSE60655.pos <- subset(GSE60655.sig, logFC > 0); nrow(GSE60655.pos)  # 16163 (p < 0.01: 3726)
GSE60655.neg <- subset(GSE60655.sig, logFC < 0); nrow(GSE60655.neg)  # 16459 (p < 0.01: 2911)

GSE114763.sig <- resultsSig[["resGSE114763p005"]]
GSE114763.pos <- subset(GSE114763.sig, logFC > 0); nrow(GSE114763.pos)  # 16745 (p < 0.01: 2513)
GSE114763.neg <- subset(GSE114763.sig, logFC < 0); nrow(GSE114763.neg)  # 18058 (p < 0.01: 2766)


#. Extracting overlaps ---- 
overlap.cpgs.pos <- intersect(GSE60655.pos$Name, GSE114763.pos$Name); length(overlap.cpgs.pos)  # 230 (p < 0.01: 10)
overlap.cpgs.neg <- intersect(GSE60655.neg$Name, GSE114763.neg$Name); length(overlap.cpgs.neg)  # 204 (p < 0.01: 15)


#### res001 (optional) ####
# GSE60655.sig <- resultsSig[["resGSE60655p001"]]
# GSE60655.pos <- subset(GSE60655.sig, logFC > 0); nrow(GSE60655.pos)  # 16163 (p < 0.01: 3726)
# GSE60655.neg <- subset(GSE60655.sig, logFC < 0); nrow(GSE60655.neg)  # 16459 (p < 0.01: 2911)
# 
# GSE114763.sig <- resultsSig[["resGSE114763p001"]]
# GSE114763.pos <- subset(GSE114763.sig, logFC > 0); nrow(GSE114763.pos)  # 16745 (p < 0.01: 2513)
# GSE114763.neg <- subset(GSE114763.sig, logFC < 0); nrow(GSE114763.neg)  # 18058 (p < 0.01: 2766)
# 
# 
# #. Extracting overlaps ---- 
# overlap.cpgs.pos <- intersect(GSE60655.pos$Name, GSE114763.pos$Name); length(overlap.cpgs.pos)  # 230 (p < 0.01: NULL)
# overlap.cpgs.neg <- intersect(GSE60655.neg$Name, GSE114763.neg$Name); length(overlap.cpgs.neg)  # 204 (p < 0.01: NULL)
#### (skip) ####


# OverlapCpGsを解析するためにマスターコードを変更(12/26/20)
# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name
length(all)

par(mfrow = c(1,2))
gst_P <- gometh(sig.cpg = OLCpGs_P, all.cpg = all, plot.bias = TRUE)
gst_N <- gometh(sig.cpg = OLCpGs_N, all.cpg = all, plot.bias = TRUE)
# warning regarding multiple symbols will always be displayed as
# there are genes that have more than one alias; not a cause for concern.

#### Top 10 GO Categories ####
# Get all the CpG sites used in the analysis to form the background
all <- DMPs[[1]]$Name
length(all) # 385220
all2 <- DMPs[[2]]$Name
length(all2) # 776698

gst.pos <- gometh(sig.cpg = overlap.cpgs.pos, plot.bias = TRUE); nrow(gst.pos)  # 22751  # 8/25 del "all.cpg = all (default: NULL)"
gst.neg <- gometh(sig.cpg = overlap.cpgs.neg, plot.bias = TRUE); nrow(gst.neg)  # 22751
gst.pos.all <- gometh(sig.cpg = overlap.cpgs.pos, all.cpg = all, plot.bias = TRUE); nrow(gst.pos)  # 22751  # 8/25 del "all.cpg = all (default: NULL)"
gst.neg.all <- gometh(sig.cpg = overlap.cpgs.neg, all.cpg = all, plot.bias = TRUE); nrow(gst.neg)  # 22751
# warning regarding multiple symbols will always be displayed as
# there are genes that have more than one alias; not a cause for concern.

t10GO.pos <- topGSA(gst.pos.all, num = Inf)
t10GO.pos <- t10GO.pos[t10GO.pos$P.DE < 0.05,]; length(t10GO.pos$TERM)  # 136
t10GO.neg <- topGSA(gst.neg.all, num = Inf)
t10GO.neg <- t10GO.neg[t10GO.neg$P.DE < 0.05,]; length(t10GO.neg$TERM)  # 126

t10GO.pos$TERM[grep("pressure", t10GO.pos$TERM)]
t10GO.neg$TERM[grep("blood pressure|cardio vascular", t10GO.neg$TERM)]

t10GO <- rbind(t10GO.pos, t10GO.neg); nrow(t10GO)  # 262
# write.csv(t10GO, "go.p005.210825.csv")
