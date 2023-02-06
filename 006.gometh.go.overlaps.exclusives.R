#### Gene Ontology mlseting (gometh) ####
# project.dir <- "/Users/Emma/GitHub/EpigenProject_Dec2020/"
# setwd(project.dir)
# source("000.util.R")
# source("001.data.prep.R")
# source("002.pdata.sort.R")
# source("003.quality.control.R")
# source("004.normalizations.R")
# source("005.differential.methylated.analysis.R")


#. significant CpG sites at FDR < 0.05 ----
# library(missMethyl)
cpgs <- list()
cpg.names <- list()
overlap <- list()
exclusive <- list()
data.type <- c("all", "sig", "sig.pos", "sig.neg")
for (i in seq_along(gseid)) {
    cat("#####", gsedir[[i]], "####", "\n")
    gse <- gseid[[i]]
    dmps <- DMPs[[i]]
    dmps.005 <- DMPs005[[i]]
    
    #### Taking overlaps of GSE60655 & GSE114763 ####
    #. loading GSE60655 data (endurance) ----
    cat("nrow(cpgs.all)", nrow(dmps), "\n")
    cat("nrow(cpgs.sig)", nrow(dmps.005), "\n")
    cpgs.sig.pos <- subset(dmps.005, logFC > 0); cat("nrow(cpgs.sig.pos)",nrow(cpgs.sig.pos), "\n")
    cpgs.sig.neg <- subset(dmps.005, logFC < 0); cat("nrow(cpgs.sig.neg)", nrow(cpgs.sig.neg), "\n\n")
    cpgs[[gse]] <- list(dmps, dmps.005, cpgs.sig.pos, cpgs.sig.neg)
    names(cpgs[[gse]]) <- data.type
    
    #. extracting cpg names ----
    cpg.names[[gse]] <- sapply(cpgs[[gse]], "[", "Name")
    
    if (i != length(gsedir)) { next }
    
    #. taking overlaps ----
    cat("#### Overlaps ####", "\n")
    for (type in names(cpg.names[[gse]])) {
        overlap[[type]] <- intersect(cpg.names[["GSE60655"]][[type]], cpg.names[["GSE114763"]][[type]])
        cat("length:", type, length(overlap[[type]]) , "\n")
    }
    
    cat("\n")
    
    #. exclusive genes ----
    cat("#### Exclusive genes ####", "\n")
    for (type in names(cpg.names[[gse]])) {
        cat("####", type, "####", "\n")
        ex.60655 <- paste0(type, ".GSE60655")
        ex.114763 <- paste0(type, ".GSE114763")
        exclusive[[ex.60655]] <- setdiff(cpg.names[["GSE60655"]][[type]], cpg.names[["GSE114763"]][[type]])
        exclusive[[ex.114763]] <- setdiff(cpg.names[["GSE114763"]][[type]], cpg.names[["GSE60655"]][[type]])
        cat("length (GSE60655):", length(exclusive[[ex.60655]]) , "\n")
        cat("length (GSE114763):", length(exclusive[[ex.114763]]) , "\n\n")
    }
}
# ##### GSE60655_End #### 
# nrow(cpgs.all) 484534 
# nrow(cpgs.sig) 47752 
# nrow(cpgs.sig.pos) 27163 
# nrow(cpgs.sig.neg) 20589 
# 
# ##### GSE114763_Res #### 
# nrow(cpgs.all) 451120 
# nrow(cpgs.sig) 27816 
# nrow(cpgs.sig.pos) 11145 
# nrow(cpgs.sig.neg) 16671 
# 
# #### Overlaps #### 
# length: all.Name 450535 
# length: sig.Name 3817 
# length: sig.pos.Name 1519 
# length: sig.neg.Name 1584 
# 
# #### Exclusive genes #### 
# #### all.Name #### 
# length (GSE60655): 33999 
# length (GSE114763): 585 
# 
# #### sig.Name #### 
# length (GSE60655): 43935 
# length (GSE114763): 23999 
# 
# #### sig.pos.Name #### 
# length (GSE60655): 25644 
# length (GSE114763): 9626 
# 
# #### sig.neg.Name #### 
# length (GSE60655): 19005 
# length (GSE114763): 15087 



#### Enrichment with GOmeth ####
gometh <- list()
go.results <- list()
gometh.interested <- c(cpg.names[["GSE60655"]][2], cpg.names[["GSE114763"]][2], 
                       overlap[2:4], exclusive[3:4], exclusive[6], exclusive[8])
names(gometh.interested)
# [1] "sig.Name"               "sig.Name"               "sig.Name"               "sig.pos.Name"           "sig.neg.Name"          
# [6] "sig.Name.GSE60655"      "sig.Name.GSE114763"     "sig.pos.Name.GSE114763" "sig.neg.Name.GSE114763"

names(gometh.interested) <- c("GSE60655sig", "GSE114763sig", 
                              "overlap.sig", "overlap.sig.pos", "overlap.sig.neg",
                              "exclusive.sig.GSE60655", "exclusive.sig.GSE114763",
                              "exclusive.sig.pos.GSE114763", "exclusive.sig.neg.GSE114763")
for (i in seq_along(gometh.interested)) {
    name <- names(gometh.interested)[i]
    cat("#### gometh():", name, "####", "\n")
    
    sig.cpg <- gometh.interested[[i]]
    gst <- gometh(sig.cpg = sig.cpg, all.cpg = overlap[["all.Name"]], plot.bias = F)
    gometh[[name]] <- gst
    
    
    #. categories ----
    cat("extracting GO terms...", "\n")
    topgsa <- topGSA(gst, num=Inf)
    topgsa <- topgsa[topgsa$P.DE < 0.05,]
    cat("P.DE < 0.05:", length(topgsa$TERM), "\n")
    go.results[[name]] <- topgsa
    
    
    #### Search the terms of interest ####
    #. positive logFC ----
    print(topgsa$TERM[grep("blood pressure|cardio vascular|inflamation", topgsa$TERM)])
    cat("position: ", grep("blood pressure|cardio vascular|inflamation", topgsa$TERM), "\n\n")
    
    # date <- format(Sys.time(), "%m.%d.%Y")
    setwd(project.dir)
    write.csv(topgsa, paste0("go.", name, ".", date, ".csv"))
    
}



#### Visualization ####
for (i in seq_along(go.results)) {
    cat("#### gometh results:", names(go.results[i]), "####", "\n")
    plot.df <- go.results[[i]] %>% 
        mutate(nlog10P = -log10(P.DE)) %>% 
        group_by(ONTOLOGY) %>% 
        slice_min(order_by = P.DE, n = 10) %>% 
        as.data.frame()
    
    print(head(plot.df, 10)); cat("\n")
    
    #. barplot ----
    idx <- order(plot.df$nlog10P, decreasing = TRUE)
    plot.df$ONTOLOGY <- factor(plot.df$ONTOLOGY)
    plot.df$TERM <- factor(plot.df$TERM, levels=rev(unique(plot.df$TERM[idx])))
    
    print(
        ggplot(plot.df, aes_string(x = "nlog10P", y = "TERM", fill = "nlog10P")) +
            theme_dose(font.size = 12) +
            scale_fill_continuous(low="blue", high="red", name = bquote(~-Log[10]~italic(P)))+
            # guide=guide_colorbar(reverse=TRUE))+
            geom_col() + # geom_bar(stat = "identity") + coord_flip() +
            # scale_y_discrete(labels = label_func) +
            facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")+
            ggtitle(paste0("GOmeth results ", names(go.results[i]))) +
            xlab(bquote(~-Log[10]~italic(P))) + ylab(NULL) +
            theme(panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank())
    )
}


cat("\n\n")

