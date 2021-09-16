#### Util variables ####
setwd("/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/")
gseid <- c("GSE60655", "GSE114763")
gsedir <- c("GSE60655_End", "GSE114763_Res")
project.dir <- "/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/"
gsedir.path <- c(paste0(project.dir, "GSE60655_End"), paste0(project.dir, "GSE114763_Res"))



#### Libraries ####
suppressMessages(
    suppressWarnings({
        library(GEOquery)
        library(lumi)  # MethyLumi
        library(wateRmelon)  # BMIQ
        library(sva)  # ComBat
        library(limma)
        library(missMethyl)
        library(RColorBrewer)
        library(stringr)
        library(dplyr)
        library(forcats)
        library(EnhancedVolcano)
        library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
        library(IlluminaHumanMethylation450kmanifest)
        library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
        library(IlluminaHumanMethylationEPICmanifest)
        library(DOSE)
    })
)

pal <- brewer.pal(8,"Dark2")  # for plotting


cat("\n\n")