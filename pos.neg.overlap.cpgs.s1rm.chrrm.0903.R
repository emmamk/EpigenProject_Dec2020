#### Gene Ontology mlseting (gometh) ####
# final test with chrremain, male files: 0903
setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/test/0902')

#### Taking overlaps of GSE60655 & GSE114763 ####
#. loading GSE60655 data (endurance) ----
end.all <- read.csv("GSE60655_End/DMPs.chrrm.lumi.clradj.bmiq.combat.0903.csv", header = TRUE); nrow(end.all)  # 472932
end.sig <- subset(end.all, P.Value < 0.05); nrow(end.sig)  # 46967
end.sig.pos <- subset(end.sig, logFC > 0); nrow(end.sig.pos)  # 20173
end.sig.pos.names <- end.sig.pos$Name
end.sig.neg <- subset(end.sig, logFC < 0); nrow(end.sig.neg)  # 26794
end.sig.neg.names <- end.sig.neg$Name

#. loading GSE114763 data (resistance) ----
res.all <- read.csv("GSE114763_Res/DMPs.s1rm.chrrm.lumi.clradj.bmiq.combat.0903.csv", header = TRUE); nrow(res.all)  # 440973
res.sig <- subset(res.all, P.Value < 0.05); nrow(res.sig)  # 27319
res.sig.pos <- subset(res.sig, logFC > 0); nrow(res.sig.pos)  # 16407
res.sig.pos.names <- res.sig.pos$Name
res.sig.neg <- subset(res.sig, logFC < 0); nrow(res.sig.neg)  # 10912
res.sig.neg.names <- res.sig.neg$Name

#. overlaps ----
overlap.cpgs.pos.p005.names <- intersect(end.sig.pos.names, res.sig.pos.names); length(overlap.cpgs.pos.p005.names)  # 1549
overlap.cpgs.neg.p005.names <- intersect(end.sig.neg.names, res.sig.neg.names); length(overlap.cpgs.neg.p005.names)  # 1513
all.cpgs.names <- intersect(end.all$Name, res.all$Name); length(all.cpgs.names)  # 440413

#. gometh() ----
suppressMessages(suppressWarnings({
    library(missMethyl)}))
gst.pos <- gometh(sig.cpg=overlap.cpgs.pos.p005.names, all.cpg=all.cpgs.names, plot.bias=F)
gst.neg <- gometh(sig.cpg=overlap.cpgs.neg.p005.names, all.cpg=all.cpgs.names, plot.bias=F)

#. categories ----
t10GO.pos <- topGSA(gst.pos, num=Inf)
t10GO.pos <- t10GO.pos[t10GO.pos$P.DE < 0.05,]; length(t10GO.pos$TERM)  # 584
t10GO.neg <- topGSA(gst.neg, num=Inf)
t10GO.neg <- t10GO.neg[t10GO.neg$P.DE < 0.05,]; length(t10GO.neg$TERM)  # 540



#### Search the terms of interest ####
#. positive logFC ----
t10GO.pos$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)  # 0

#. negative logFC ----
t10GO.neg$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.neg$TERM)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO.neg$TERM)  # 0


write.csv(t10GO.pos, "overlap.go.p005.pos.s1rm.chrrm.lumi.clradj.bmiq.combat.0903.csv")
write.csv(t10GO.neg, "overlap.go.p005.neg.s1rm.chrrm.lumi.clradj.bmiq.combat.0903.csv")
t10GO.pos$TERM[1:10]
#  [1] "adaptive immune response"                    "immune system development"                   "sensory organ development"                  
# [4] "sensory organ morphogenesis"                 "hematopoietic or lymphoid organ development" "cell fate commitment"                       
# [7] "cell-cell adhesion"                          "alpha-beta T cell differentiation"           "cell adhesion"                              
# [10] "biological adhesion"       

t10GO.neg$TERM[1:10]
# [1] "myofibril"                                     "contractile fiber"                             "sarcomere"                                    
# [4] "muscle system process"                         "I band"                                        "A band"                                       
# [7] "muscle contraction"                            "negative regulation of protein polymerization" "Z disc"                                       
# [10] "myofibril assembly" 

