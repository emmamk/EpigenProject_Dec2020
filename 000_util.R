#### Util variables ####
setwd("/Users/Emma/Documents/Bioinformatics/Epigenetics/")
gseid <- c("GSE60655", "GSE114763")
gsedir <- c("GSE60655_End", "GSE114763_Res")

#### Libraries ####
suppressMessages({
    library(GEOquery)
    library(minfi)
    library(stringr)
    library(dplyr)
    library(readxl)
    library(DMRcate)
    library(limma)
    library(RColorBrewer)
    library(missMethyl)
    library(minfiData)
    library(Gviz)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
    library(IlluminaHumanMethylationEPICmanifest)
    })
