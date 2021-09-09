#### Gene Ontology mlseting (gometh) ####
# final test with chrremain, male files: 0903
setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/')

# loading libraries ----
suppressMessages(
    suppressWarnings({
        library(missMethyl)
        library(enrichplot)
        library(ggplot2)
        library(dplyr)
        library(DOSE)
        # library(clusterProfiler)
    })
    )


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
res.sig.pos <- subset(res.sig, logFC > 0); nrow(res.sig.pos)  # 11152
res.sig.neg <- subset(res.sig, logFC < 0); nrow(res.sig.neg)  # 16678
res.sig.names <- res.sig$Name
res.sig.pos.names <- res.sig.pos$Name
res.sig.neg.names <- res.sig.neg$Name

#. overlaps ----
overlap.cpgs.p005.names <- intersect(end.sig.names, res.sig.names); length(overlap.cpgs.p005.names)  # 3817
overlap.cpgs.pos.p005.names <- intersect(end.sig.pos.names, res.sig.pos.names); length(overlap.cpgs.pos.p005.names)  # 1520
overlap.cpgs.neg.p005.names <- intersect(end.sig.neg.names, res.sig.neg.names); length(overlap.cpgs.neg.p005.names)  # 1584
all.cpgs.names <- intersect(end.all$Name, res.all$Name); length(all.cpgs.names)  # 450906

#. gometh() ----
gst <- missMethyl::gometh(sig.cpg=overlap.cpgs.p005.names, all.cpg=all.cpgs.names, plot.bias=F)
gst.pos <- gometh(sig.cpg=overlap.cpgs.pos.p005.names, all.cpg=all.cpgs.names, plot.bias=F)
gst.neg <- gometh(sig.cpg=overlap.cpgs.neg.p005.names, all.cpg=all.cpgs.names, plot.bias=F)


#. categories ----
t10GO <- topGSA(gst, num=Inf)
t10GO <- t10GO[t10GO$P.DE < 0.05,]; length(t10GO$TERM); head(t10GO)  # 407
#            ONTOLOGY                         TERM   N  DE         P.DE FDR
# GO:0048736       BP        appendage development 180  51 0.0002972569   1
# GO:0060173       BP             limb development 180  51 0.0002972569   1
# GO:1990845       BP       adaptive thermogenesis 155  40 0.0003042782   1
# GO:0030029       BP actin filament-based process 797 160 0.0003380761   1
# GO:0030016       CC                    myofibril 225  52 0.0003613794   1
# GO:0043292       CC            contractile fiber 236  53 0.0004609175   1


t10GO.pos <- topGSA(gst.pos, num=Inf)
t10GO.pos <- t10GO.pos[t10GO.pos$P.DE < 0.05,]; length(t10GO.pos$TERM); head(t10GO.pos)  # 536
#            ONTOLOGY                  TERM   N DE         P.DE          FDR
# GO:0030016       CC             myofibril 225 40 1.290130e-08 0.0002935174
# GO:0043292       CC     contractile fiber 236 40 4.248486e-08 0.0004832865
# GO:0030017       CC             sarcomere 205 36 1.017755e-07 0.0007718318
# GO:0003012       BP muscle system process 464 58 2.042246e-06 0.0116157822
# GO:0031674       CC                I band 139 26 8.198406e-06 0.0373043883
# GO:0031672       CC                A band  37 11 2.400576e-05 0.0910258564


t10GO.neg <- topGSA(gst.neg, num=Inf)
t10GO.neg <- t10GO.neg[t10GO.neg$P.DE < 0.05,]; length(t10GO.neg$TERM); head(t10GO.neg)  # 597
#            ONTOLOGY                                        TERM    N DE         P.DE       FDR
# GO:0002250       BP                    adaptive immune response  403 39 3.007276e-05 0.6193752
# GO:0002520       BP                   immune system development 1003 84 1.024987e-04 0.6193752
# GO:0007423       BP                   sensory organ development  574 62 1.118460e-04 0.6193752
# GO:0090596       BP                 sensory organ morphogenesis  268 36 1.243848e-04 0.6193752
# GO:0048534       BP hematopoietic or lymphoid organ development  949 80 1.631535e-04 0.6193752
# GO:0045165       BP                        cell fate commitment  278 38 1.766110e-04 0.6193752



#### Search the terms of interest ####
#. all ----
t10GO$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO$TERM)]  # 0
grep("blood pressure|cardio vascular|inflamation", t10GO$TERM)  # 0

#. positive logFC ----
t10GO.pos$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)]
# [1] "nervous system process involved in regulation of systemic arterial blood pressure"
grep("blood pressure|cardio vascular|inflamation", t10GO.pos$TERM)  # 465

#. negative logFC ----
t10GO.neg$TERM[grep("blood pressure|cardio vascular|inflamation", t10GO.neg$TERM)]
grep("blood pressure|cardio vascular|inflamation", t10GO.neg$TERM)  # 0


# write.csv(t10GO, "overlap.go.p005.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")
# write.csv(t10GO.pos, "overlap.go.p005.pos.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")
# write.csv(t10GO.neg, "overlap.go.p005.neg.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv")



#### Visualization ####
plot.df <- t10GO.pos %>% 
    mutate(GeneRatio = DE / N, nlog10P = -log10(P.DE)) %>% 
    group_by(ONTOLOGY) %>% 
    slice_min(order_by = P.DE, n = 10) %>% 
    as.data.frame()
head(plot.df)
summary(plot.df$GeneRatio)

#. dotplot ----
idx <- order(plot.df$GeneRatio, decreasing = TRUE)
plot.df$ONTOLOGY <- factor(plot.df$ONTOLOGY)
plot.df$TERM <- factor(plot.df$TERM, levels=rev(unique(plot.df$TERM[idx])))

ggplot(plot.df, aes_string(plot.df$GeneRatio, y = plot.df$TERM, size=plot.df$DE, color=plot.df$P.DE)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = "p-value",
                           guide=guide_colorbar(reverse=TRUE)) +
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")+
    # scale_y_discrete(labels = test$TERM)+
    ylab(NULL) + xlab("GeneRatio") + ggtitle("title") + theme_dose(font.size = 12) +
    scale_size(range=c(3, 8))


#. barplot ----
idx <- order(plot.df$logP, decreasing = TRUE)
plot.df$ONTOLOGY <- factor(plot.df$ONTOLOGY)
plot.df$TERM <- factor(plot.df$TERM, levels=rev(unique(plot.df$TERM[idx])))

ggplot(plot.df, aes_string(x = "logP", y = "TERM", fill = "logP")) +
    theme_dose(font.size = 12) +
    scale_fill_continuous(low="blue", high="red", name = bquote(~-Log[10]~italic(P)))+
                          # guide=guide_colorbar(reverse=TRUE))+
    geom_col() + # geom_bar(stat = "identity") + coord_flip() +
    # scale_y_discrete(labels = label_func) +
    facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")+
    ggtitle("title") + xlab(bquote(~-Log[10]~italic(P))) + ylab(NULL)

