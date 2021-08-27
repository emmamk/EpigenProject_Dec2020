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

# names(resultsSig)
# "resGSE60655adj02"  "resGSE114763adj02" "resGSE60655p005"   "resGSE114763p005" 
# "resGSE60655p001"  "resGSE114763p001"

GSE60655.sig.p005 <- resultsSig[["resGSE60655p005"]]
# GSE60655.pos <- subset(GSE60655.sig.p005, logFC > 0); nrow(GSE60655.pos)  # 19875
# GSE60655.neg <- subset(GSE60655.sig.p005, logFC < 0); nrow(GSE60655.neg)  # 20168

GSE114763.sig.p005 <- resultsSig[["resGSE114763p005"]]
# GSE114763.pos <- subset(GSE114763.sig.p005, logFC > 0); nrow(GSE114763.pos)  # 21640
# GSE114763.neg <- subset(GSE114763.sig.p005, logFC < 0); nrow(GSE114763.neg)  # 21311


#. Extracting overlaps ---- 
# merging with CpG site names
innerjoin005 <- inner_join(GSE60655.sig.p005, GSE114763.sig.p005, by = "Name"); nrow(innerjoin005)  # 1821 (2281)
ij005names <- innerjoin005$Name
is005names <- intersect(GSE60655.sig.p005$Name, GSE114763.sig.p005$Name); length(is005names)  # 1821 (2281)

names(DMPs)  # "GSE60655"  "GSE114763"
ij.all.cpgs <- inner_join(DMPs[["GSE60655"]], DMPs[["GSE114763"]], by = "Name"); nrow(ij.all.cpgs)  # 357086 (434778)
ij.all.names <- ij.all.cpgs$Name
is.all.names <- intersect(DMPs[["GSE60655"]]$Name, DMPs[["GSE114763"]]$Name); length(is.all.names)  # 357086 (434778)

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



#### Top 10 GO Categories ####
# Get all the CpG sites used in the analysis to form the background
# all <- DMPs[[1]]$Name

# all <- all.cpgs$Name
gst.ij005 <- gometh(sig.cpg = ij005names, all.cpg = ij.all.names, plot.bias = F); nrow(gst.ij005)  # 22732 (22751)  # memo: same results wiyh array.type = "EPIC"
gst.is005 <- gometh(sig.cpg = is005names, all.cpg = is.all.names, plot.bias = F); nrow(gst.is005)  # 22732
# gst.pos <- gometh(sig.cpg = overlap.cpgs.pos, plot.bias = TRUE); nrow(gst.pos)  # 22751  # 8/25 del "all.cpg = all (default: NULL)"
# gst.neg <- gometh(sig.cpg = overlap.cpgs.neg, plot.bias = TRUE); nrow(gst.neg)  # 22751
# gst.pos.all <- gometh(sig.cpg = overlap.cpgs.pos, all.cpg = all, plot.bias = TRUE); nrow(gst.pos)  # 22751  # 8/25 del "all.cpg = all (default: NULL)"
# gst.neg.all <- gometh(sig.cpg = overlap.cpgs.neg, all.cpg = all, plot.bias = TRUE); nrow(gst.neg)  # 22751
# warning regarding multiple symbols will always be displayed as
# there are genes that have more than one alias; not a cause for concern.

# t10GO.pos <- topGSA(gst.pos.all, num = Inf)
# t10GO.pos <- t10GO.pos[t10GO.pos$P.DE < 0.05,]; length(t10GO.pos$TERM)  # 136
# t10GO.neg <- topGSA(gst.neg.all, num = Inf)
# t10GO.neg <- t10GO.neg[t10GO.neg$P.DE < 0.05,]; length(t10GO.neg$TERM)  # 126
ij.t10GO.005 <- topGSA(gst.ij005, num = Inf)
ij.t10GO.005 <- ij.t10GO.005[ij.t10GO.005$P.DE < 0.05,]; length(ij.t10GO.005$TERM)  # 352
ij.t10GO.005$TERM[1:20]

is.t10GO.005 <- topGSA(gst.is005, num = Inf)
is.t10GO.005 <- is.t10GO.005[is.t10GO.005$P.DE < 0.05,]; length(is.t10GO.005$TERM)  # 352
is.t10GO.005$TERM[1:20]

ij.t10GO.005$TERM[grep("blood pressure|cardio vascular|inflamation", ij.t10GO.005$TERM)]
# [1] "angiotensin-mediated vasodilation involved in regulation of systemic arterial blood pressure"
# ij.t10GO.005[grep("blood pressure|cardio vascular|inflamation", ij.t10GO.005$TERM),]
grep("blood pressure|cardio vascular|inflamation", ij.t10GO.005$TERM)  # 77
    
is.t10GO.005$TERM[grep("blood pressure|cardio vascular|inflamation", is.t10GO.005$TERM)]
grep("blood pressure|cardio vascular|inflamation", is.t10GO.005$TERM)  # 77
    
# "angiotensin-mediated vasodilation involved in regulation of systemic arterial blood pressure"

# t10GO.pos$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)]
# t10GO.neg$TERM[grep("blood pressure|cardio vascular", t10GO.neg$TERM)]

# t10GO <- rbind(t10GO.pos, t10GO.neg); nrow(t10GO)  # 262
# write.csv(t10GO, "go.p005.210825.csv")



### memo: condition -- no removal of chrX/Y, subsetted with male subjects
### ij.t10GO.005$TERM[1:20]
# [1] "error-free translesion synthesis"                       
# [2] "malate dehydrogenase (decarboxylating) (NAD+) activity" 
# [3] "malate dehydrogenase (decarboxylating) (NADP+) activity"
# [4] "sex determination"                                      
# [5] "MLL3/4 complex"                                         
# [6] "oxaloacetate decarboxylase activity"                    
# [7] "malic enzyme activity"                                  
# [8] "rRNA binding"                                           
# [9] "glossopharyngeal nerve development"                     
# [10] "regulation of mitotic nuclear division"                 
# [11] "maintenance of sister chromatid cohesion"               
# [12] "maintenance of mitotic sister chromatid cohesion"       
# [13] "positive regulation of metanephros development"         
# [14] "PML body"                                               
# [15] "regulation of cell adhesion"                            
# [16] "vagus nerve development"                                
# [17] "establishment of protein localization to chromosome"    
# [18] "positive regulation of cell adhesion"                   
# [19] "positive regulation of protein kinase B signaling"      
# [20] "germ cell migration" 
