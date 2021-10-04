#### Gene Ontology by gometh(KEGG) ####
# setwd("/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/")
# source("000.util.R")
# source("001.data.prep.R")
# source("002.pdata.sort.R")
# source("003.quality.control.R")
# source("004.normalizations.R")
# source("005.differential.methylated.analysis.R")
# source(006.gometh.go.overlaps.exclusives.R)

#### Enrichment with GOmeth (KEGG) ####
kegg <- list()
kegg.results <- list()
kegg.interested <- c(cpg.names[["GSE60655"]][2], cpg.names[["GSE114763"]][2],
                     overlap[2], exclusive[3:4])
names(kegg.interested)  # [1] "sig.Name"  "sig.Name"  "sig.Name"  "sig.Name.GSE60655"  "sig.Name.GSE114763"
names(kegg.interested) <- c("GSE60655sig", "GSE114763sig", "overlap.sig", "exclusive.sig.GSE60655", "exclusive.sig.GSE114763")
for (i in seq_along(kegg.interested)) {
    name <- names(kegg.interested)[i]
    cat("#### gometh(kegg):", name, "####", "\n")
    
    sig.cpg <- kegg.interested[[i]]
    gst <- gometh(sig.cpg = sig.cpg, all.cpg = overlap[["all.Name"]], plot.bias = F, collection = "KEGG")
    kegg[[name]] <- gst
    
    
    #. categories ----
    cat("extracting KEGG terms...", "\n")
    topgsa <- topGSA(gst, num=Inf)
    topgsa <- topgsa[topgsa$P.DE < 0.05,]
    cat("P.DE < 0.05:", length(topgsa$Description), "\n")
    kegg.results[[name]] <- topgsa
    
    
    #### Search the terms of interest ####
    #. positive logFC ----
    print(topgsa$Description[grep("blood pressure|cardio vascular|inflamation", topgsa$Description)])
    cat("position: ", grep("blood pressure|cardio vascular|inflamation", topgsa$Description), "\n\n")
    
    date <- format(Sys.time(), "%m.%d.%Y")
    setwd(project.dir)
    write.csv(topgsa, paste0("kegg.", name, ".", date, ".csv"))
    
}



#### Visualization ####
for (i in seq_along(kegg.results)) {
    res.name = names(kegg.results[i])
    cat("#### dotplot for kegg:", res.name, "####", "\n")
    
    #. dotplot ----
    if (res.name == "exclusive.sig.GSE60655") {
        for (j in c(20, 30, 40)) {
            plot.df <- kegg.results[[i]]
            plot.df <- plot.df %>% 
                slice_max(order_by = DE, n = j)
            idx <- order(plot.df$DE, decreasing = TRUE)
            plot.df$Description <- factor(plot.df$Description, levels=rev(unique(plot.df$Description[idx])))
            
            png(paste0(plot.dir, "kegg.dotplot.", res.name, ".", j, ".", date, ".png"),
                width = 1280, height = 1000)
            print(
            ggplot(plot.df, 
                   aes_string(plot.df$DE, y = plot.df$Description, size=plot.df$DE, color=plot.df$P.DE)) +
                geom_point() +
                scale_color_continuous(low="red", high="blue", name = "p-value",
                                       guide=guide_colorbar(reverse=TRUE)) +
                # facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")+
                # scale_y_discrete(labels = test$TERM)+
                ylab(NULL) + xlab("Count") + ggtitle(paste0("KEGG ", res.name, " ", j)) +
                theme_dose(font.size = 16) +
                scale_size(range=c(3, 8))
            )
            dev.off()
        }
    } else {
        plot.df <- kegg.results[[i]]
        idx <- order(plot.df$DE, decreasing = TRUE)
        plot.df$Description <- factor(plot.df$Description, levels=rev(unique(plot.df$Description[idx])))
        
        png(paste0(plot.dir, "kegg.dotplot.", res.name, ".", date, ".png"), width = 1280, height = 1000)
        print(
            ggplot(plot.df, 
                   aes_string(plot.df$DE, y = plot.df$Description, size=plot.df$DE, color=plot.df$P.DE)) +
                geom_point() +
                scale_color_continuous(low="red", high="blue", name = "p-value",
                                       guide=guide_colorbar(reverse=TRUE)) +
                # facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")+
                # scale_y_discrete(labels = test$TERM)+
                ylab(NULL) + xlab("Count") + ggtitle(paste0("KEGG ", res.name)) +
                theme_dose(font.size = 16) +
                scale_size(range=c(3, 8))
        )
        dev.off()
    }

        
    #. barplot ----
    # plot.df <- kegg.results[[i]] %>%
    #     mutate(nlog10P = -log10(P.DE)) %>%
    #     group_by(ONTOLOGY) %>%
    #     slice_min(order_by = P.DE, n = 10) %>%
    #     as.data.frame()
    # 
    # print(head(plot.df, 10)); cat("\n")
    
    # idx <- order(plot.df$nlog10P, decreasing = TRUE)
    # plot.df$ONTOLOGY <- factor(plot.df$ONTOLOGY)
    # plot.df$TERM <- factor(plot.df$TERM, levels=rev(unique(plot.df$TERM[idx])))
    # 
    # print(
    #     ggplot(plot.df, aes_string(x = "nlog10P", y = "TERM", fill = "nlog10P")) +
    #         theme_dose(font.size = 12) +
    #         scale_fill_continuous(low="blue", high="red", name = bquote(~-Log[10]~italic(P)))+
    #         # guide=guide_colorbar(reverse=TRUE))+
    #         geom_col() + # geom_bar(stat = "identity") + coord_flip() +
    #         # scale_y_discrete(labels = label_func) +
    #         facet_grid(ONTOLOGY ~ ., scales = "free", space = "free")+
    #         ggtitle(paste0("GOmeth results ", names(go.results[i]))) +
    #         xlab(bquote(~-Log[10]~italic(P))) + ylab(NULL) +
    #         theme(panel.grid.major = element_blank(),
    #               panel.grid.minor = element_blank())
    # )
}

cat("\n\n")
