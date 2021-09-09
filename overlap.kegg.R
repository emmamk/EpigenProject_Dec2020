#### Gene Ontology mlseting (gometh) ####
# final test with chrremain, male files: 0903
setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/test/0902ref.before')

#### Taking overlaps of GSE60655 & GSE114763 ####
#. loading GSE60655 data (endurance) ----
end.all <- read.csv("GSE60655_End/DMPs.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv", header = TRUE); nrow(end.all)  # 484534
end.sig <- subset(end.all, P.Value < 0.05); nrow(end.sig)  # 47752
end.sig.names <- end.sig$Name
end.sig.pos <- subset(end.sig, logFC > 0); nrow(end.sig.pos)  # 27163
end.sig.pos.names <- end.sig.pos$Name
end.sig.neg <- subset(end.sig, logFC < 0); nrow(end.sig.neg)  # 20589
end.sig.neg.names <- end.sig.neg$Name

#. loading GSE114763 data (resistance) ----
res.all <- read.csv("GSE114763_Res/DMPs.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv", header = TRUE); nrow(res.all)  # 451492
res.sig <- subset(res.all, P.Value < 0.05); nrow(res.sig)  # 27830
res.sig.names <- res.sig$Name
res.sig.pos <- subset(res.sig, logFC > 0); nrow(res.sig.pos)  # 11152
res.sig.pos.names <- res.sig.pos$Name
res.sig.neg <- subset(res.sig, logFC < 0); nrow(res.sig.neg)  # 16678
res.sig.neg.names <- res.sig.neg$Name

#. overlaps ----
overlap.cpgs.p005.names <- intersect(end.sig.names, res.sig.names); length(overlap.cpgs.p005.names)  # 3817
overlap.cpgs.pos.p005.names <- intersect(end.sig.pos.names, res.sig.pos.names); length(overlap.cpgs.pos.p005.names)  # 1520
overlap.cpgs.neg.p005.names <- intersect(end.sig.neg.names, res.sig.neg.names); length(overlap.cpgs.neg.p005.names)  # 1584
all.cpgs.names <- intersect(end.all$Name, res.all$Name); length(all.cpgs.names)  # 450906

#. gometh() ----
suppressMessages(suppressWarnings({
    library(missMethyl)}))
gst <- gometh(sig.cpg=overlap.cpgs.p005.names, all.cpg=all.cpgs.names, collection = "KEGG", plot.bias=F)
gst.pos <- gometh(sig.cpg=overlap.cpgs.pos.p005.names, all.cpg=all.cpgs.names, collection = "KEGG", plot.bias=F)
gst.neg <- gometh(sig.cpg=overlap.cpgs.neg.p005.names, all.cpg=all.cpgs.names, collection = "KEGG", plot.bias=F)

#. categories ----
t10GO <- topGSA(gst, num=Inf)
t10GO <- t10GO[t10GO$P.DE < 0.05,]; length(t10GO$Description)  # 19
t10GO.pos <- topGSA(gst.pos, num=Inf)
t10GO.pos <- t10GO.pos[t10GO.pos$P.DE < 0.05,]; length(t10GO.pos$Description)  # 14
t10GO.neg <- topGSA(gst.neg, num=Inf)
t10GO.neg <- t10GO.neg[t10GO.neg$P.DE < 0.05,]; length(t10GO.neg$Description)  # 16


#### Search the terms of interest ####
#. all ----
t10GO$Description[grep("blood pressure|cardio vascular|inflamation", t10GO$Description)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO$Description)  # 0

#. positive logFC ----
t10GO.pos$Description[grep("blood pressure|cardio vascular|inflamation", t10GO.pos$Description)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO.pos$Description)  # 0

#. negative logFC ----
t10GO.neg$Description[grep("blood pressure|cardio vascular|inflamation", t10GO.neg$Description)]
grep("blood pressure|cardio vascular|inflamation", t10GO.neg$Description)  # 0


write.csv(t10GO, "overlap.kegg.p005.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")
write.csv(t10GO.pos, "overlap.kegg.p005.pos.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")
write.csv(t10GO.neg, "overlap.kegg.p005.neg.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")


