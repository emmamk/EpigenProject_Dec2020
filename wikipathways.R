#### Gene Ontology with rWikiPathways ####
# chrremain, male, s1(GSE114763) & s27(GSE60655) removed: 0907
setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/')

suppressMessages(
    suppressWarnings({
    library("DOSE")
    library("GO.db")
    library("GSEABase")
    library("org.Hs.eg.db")
    library("clusterProfiler")
    library("dplyr")
    library("tidyr")
    library("ggplot2")
    library("stringr")
    library("RColorBrewer")
    library("rWikiPathways")
    library("enrichplot")
    library("RCy3")
    library("AnnotationDbi")
    library("minfi")
    }))


#### Taking overlaps of GSE60655 & GSE114763 ####
#. loading GSE60655 data (endurance) ----
end.all <- read.csv("GSE60655_End/DMPs.chrremain.refbefore.lumi.clradj.bmiq.combat.0904.csv", header = TRUE); nrow(end.all)  # 484534
end.sig <- subset(end.all, P.Value < 0.05)
end.sig.up <- end.sig[end.sig$logFC > 0, "UCSC_RefGene_Name"]; length(end.sig.up)  # 27163
end.sig.dn <- end.sig[end.sig$logFC < 0, "UCSC_RefGene_Name"]; length(end.sig.dn)  # 20589

#. loading GSE114763 data (resistance) ----
res.all <- read.csv("GSE114763_Res/DMPs.s1rm.chrremain.lumi.clradj.bmiq.combat.refbefore.0906.csv", header = TRUE); nrow(res.all)  # 451492
res.sig <- subset(res.all, P.Value < 0.05)
res.sig.up <- res.sig[res.sig$logFC > 0, "UCSC_RefGene_Name"]; length(res.sig.up)  # 27163
res.sig.dn <- res.sig[res.sig$logFC < 0, "UCSC_RefGene_Name"]; length(res.sig.dn)  # 20589

#. overlaps ----
overlap.all <- intersect(end.all$UCSC_RefGene_Name, res.all$UCSC_RefGene_Name); length(overlap.all)  # 41652
overlap.sig <- intersect(end.sig$UCSC_RefGene_Name, res.sig$UCSC_RefGene_Name); length(overlap.sig)  # 41652
overlap.up <- intersect(end.sig.up, res.sig.up); length(overlap.up) # 3253
overlap.dn <- intersect(end.sig.dn, res.sig.dn); length(overlap.dn) # 3309
overlap.list <- list(overlap.all, overlap.sig, overlap.up, overlap.dn)
names(overlap.list) <- c("overlap.all", "overlap.sig", "overlap.up", "overlap.dn")


#### Adding Annotation: ENSEMBL IDs ####
# library(org.Hs.eg.db)
for (i in seq_along(overlap.list)) {
    x <- overlap.list[[i]]
    anno <- AnnotationDbi::select(org.Hs.eg.db,
                          keys = x,
                          keytype = "SYMBOL",
                          column = "ENSEMBL")
    x <- data.frame(SYMBOL = x)
    x <- left_join(x, anno)
    cat("nrow:", nrow(x), "\n\n")
    
    overlap.list[[i]] <- x
    }
bkgd.genes <- overlap.list[[1]]$ENSEMBL; head(bkgd.genes)
sig.genes <- overlap.list[[2]]$ENSEMBL; head(sig.genes)
up.genes <- overlap.list[[3]]$ENSEMBL; head(up.genes)
dn.genes <- overlap.list[[4]]$ENSEMBL; head(dn.genes)

#. getting ENTREZID for future use ----
bkgd.genes.entrez <- clusterProfiler::bitr(bkgd.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db); head(bkgd.genes.entrez)
sig.genes.entrez <- bitr(sig.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db); head(sig.genes.entrez)
up.genes.entrez <- bitr(up.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db); head(up.genes.entrez)
dn.genes.entrez <- bitr(dn.genes,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db); head(dn.genes.entrez)



#### Preparing rWikiPathways terms ####
wp.hs.gmt <- rWikiPathways::downloadPathwayArchive(organism="Homo sapiens", format = "gmt")
wp2gene <- readPathwayGMT(wp.hs.gmt)
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene); head(wpid2gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name); head(wpid2name) #TERM2NAME



#### WikiPathways: enricher() ####
#. sig.genes ----
ewp.sig <- clusterProfiler::enricher(
    sig.genes.entrez[[2]],
    universe = bkgd.genes.entrez[[2]],
    pAdjustMethod = "none",
    pvalueCutoff = 0.05,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
ewp.sig <- DOSE::setReadable(ewp.sig, org.Hs.eg.db, keyType = "ENTREZID"); nrow(ewp.sig)  # p < 0.05: 58, FDR < 0.05: 22
ewp.sig$Description
# [1] "Osteoblast differentiation and related diseases"                                    "Neural crest differentiation"                                                      
# [3] "Embryonic stem cell pluripotency pathways"                                          "Wnt signaling pathway and pluripotency"                                            
# [5] "Breast cancer pathway"                                                              "Insulin signaling"                                                                 
# [7] "Hair follicle development: organogenesis - part 2 of 3"                             "Nephrogenesis"                                                                     
# [9] "Wnt signaling in kidney disease"                                                    "Myometrial relaxation and contraction pathways"                                    
# [11] "Focal adhesion: PI3K-Akt-mTOR-signaling pathway"                                    "Focal adhesion"                                                                    
# [13] "Endochondral ossification"                                                          "Endochondral ossification with skeletal dysplasias"                                
# [15] "Gastrin signaling pathway"                                                          "Canonical and non-canonical Notch signaling"                                       
# [17] "B cell receptor signaling pathway"                                                  "Hippo signaling regulation pathways"                                               
# [19] "Malignant pleural mesothelioma"                                                     "Epithelial to mesenchymal transition in colorectal cancer"                         
# [21] "4-hydroxytamoxifen, Dexamethasone, and Retinoic Acids Regulation of p27 Expression" "Ebola virus pathway in host"                                                       
# [23] "Notch signaling pathway (Netpath)"                                                  "Alzheimer's disease"                                                               
# [25] "Mesodermal commitment pathway"                                                      "Deregulation of Rab and Rab effector genes in bladder cancer"                      
# [27] "Endoderm differentiation"                                                           "Hippo-Merlin signaling dysregulation"                                              
# [29] "Heart development"                                                                  "Hedgehog signaling pathway Netpath"                                                
# [31] "Splicing factor NOVA regulated synaptic proteins"                                   "Differentiation of white and brown adipocyte"                                      
# [33] "Notch signaling"                                                                    "Ectoderm differentiation"                                                          
# [35] "Interleukin-1 (IL-1) structural pathway"                                            "PI3K-Akt signaling pathway"                                                        
# [37] "Primary focal segmental glomerulosclerosis (FSGS)"                                  "ncRNAs involved in Wnt signaling in hepatocellular carcinoma"                      
# [39] "Adipogenesis"                                                                       "Nephrotic syndrome"                                                                
# [41] "Wnt signaling"                                                                      "Non-genomic actions of 1,25 dihydroxyvitamin D3"                                   
# [43] "lncRNA in canonical Wnt signaling and colorectal cancer"                            "MAPK signaling pathway"                                                            
# [45] "MicroRNAs in cardiomyocyte hypertrophy"                                             "Integrin-mediated cell adhesion"                                                   
# [47] "Host-pathogen interaction of human coronaviruses - MAPK signaling"                  "FGFR3 signaling in chondrocyte proliferation and terminal differentiation"         
# [49] "Mammalian disorder of sexual development"                                           "Glioblastoma signaling pathways"                                                   
# [51] "Cardiac hypertrophic response"                                                      "Vitamin D in inflammatory diseases"                                                
# [53] "Small cell lung cancer"                                                             "miR-509-3p alteration of YAP1/ECM axis"                                            
# [55] "Resistin as a regulator of inflammation"                                            "TGF-beta signaling pathway"                                                        
# [57] "Angiopoietin-like protein 8 regulatory pathway"                                     "Fragile X syndrome" 

barplot(ewp.sig, showCategory = 20)
dotplot(ewp.sig, showCategory = 20)
ewp.sig.pairwise <- enrichplot::pairwise_termsim(ewp.sig)
emapplot(ewp.sig.pairwise, showCategory = 20)


#. up-regulated genes ----
ewp.up <- clusterProfiler::enricher(
    up.genes.entrez[[2]],
    universe = bkgd.genes.entrez[[2]],
    pAdjustMethod = "none",
    pvalueCutoff = 0.05,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
ewp.up <- DOSE::setReadable(ewp.up, org.Hs.eg.db, keyType = "ENTREZID"); nrow(ewp.up)  # 0
head(ewp.up)  # 0
# barplot(ewp.up, showCategory = 20)
# dotplot(ewp.up, showCategory = 20)
# emapplot(ewp.up, showCategory = 20)


#. down-regulated genes ----
ewp.dn <- enricher(
    dn.genes.entrez[[2]],
    universe = bkgd.genes.entrez[[2]],
    pAdjustMethod = "none",
    pvalueCutoff = 0.05,
    TERM2GENE = wpid2gene,
    TERM2NAME = wpid2name)
ewp.dn <- setReadable(ewp.dn, org.Hs.eg.db, keyType = "ENTREZID"); nrow(ewp.dn)  # 42
ewp.dn$Description
# [1] "Osteoblast differentiation and related diseases"              "Neural crest differentiation"                                
# [3] "Epithelial to mesenchymal transition in colorectal cancer"    "Dopaminergic neurogenesis"                                   
# [5] "Wnt signaling in kidney disease"                              "Embryonic stem cell pluripotency pathways"                   
# [7] "Hair follicle development: organogenesis - part 2 of 3"       "Nephrogenesis"                                               
# [9] "Ectoderm differentiation"                                     "Breast cancer pathway"                                       
# [11] "Notch signaling pathway (Netpath)"                            "Endochondral ossification"                                   
# [13] "Endochondral ossification with skeletal dysplasias"           "Hippo signaling regulation pathways"                         
# [15] "Wnt signaling pathway and pluripotency"                       "Hippo-Merlin signaling dysregulation"                        
# [17] "Ebola virus pathway in host"                                  "Adipogenesis"                                                
# [19] "Mammalian disorder of sexual development"                     "Endoderm differentiation"                                    
# [21] "miR-509-3p alteration of YAP1/ECM axis"                       "Hedgehog signaling pathway Netpath"                          
# [23] "Malignant pleural mesothelioma"                               "Somatic sex determination"                                   
# [25] "Allograft Rejection"                                          "Mesodermal commitment pathway"                               
# [27] "Wnt signaling"                                                "Heart development"                                           
# [29] "Canonical and non-canonical Notch signaling"                  "miRNA targets in ECM and membrane receptors"                 
# [31] "Notch signaling"                                              "Focal adhesion"                                              
# [33] "Differentiation of white and brown adipocyte"                 "Cell migration and invasion through p75NTR"                  
# [35] "ncRNAs involved in Wnt signaling in hepatocellular carcinoma" "Hedgehog signaling pathway"                                  
# [37] "Insulin signaling"                                            "Focal adhesion: PI3K-Akt-mTOR-signaling pathway"             
# [39] "15q13.3 copy number variation syndrome"                       "B cell receptor signaling pathway"                           
# [41] "Splicing factor NOVA regulated synaptic proteins"             "Netrin-UNC5B signaling pathway" 


write.csv(ewp.sig, "wikipathways.sig.p005.0907.csv", row.names = FALSE)
write.csv(ewp.dn, "wikipathways.dn.p005.0907.csv", row.names = FALSE)


#. Quick Visualization (down-regulated genes) ----


barplot(ewp.dn, showCategory = 20)
dotplot(ewp.dn, showCategory = 20)

ewp.dn.pairwise <- enrichplot::pairwise_termsim(ewp.dn)
emapplot(ewp.dn.pairwise, showCategory = 20)

