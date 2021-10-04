#### Util variables ####
setwd("/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/")
project.dir <- "/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/" # getwd()
date <- format(Sys.time(), "%m.%d.%Y")
gseid <- c("GSE60655", "GSE114763")
gsedir <- c("GSE60655_End", "GSE114763_Res")
gsedir.path <- c(paste0(project.dir, "GSE60655_End"), paste0(project.dir, "GSE114763_Res"))

if (!dir.exists(paste0(project.dir, "plot"))) {
    dir.create("plot")
}
plot.dir <- paste0(project.dir, "plot/")


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
        library(ggplot2)
    })
)

pal <- brewer.pal(8,"Dark2")  # for plotting


cat("\n\n")