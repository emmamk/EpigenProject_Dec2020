#### Removing Cross Reactive Probes (450k) ####
# (multiple places in the genome)
# Chen et al. 2013: 450k
# https://www.tandfonline.com/doi/full/10.4161/epi.23470

source("000_util.R")
# source("003_quality_control.R")

#. Obtaining Reference Files for 450k ----
# https://bioc.ism.ac.jp/packages/devel/bioc/vignettes/ramwas/inst/doc/RW5a_matrix.html#probes-with-snps-and-in-cross-reactive-regions
if (!file.exists("DataDir")) {
    dir.create("DataDir")
    }
if (!file.exists("DataDir/i450k_ns.xlsx")) {
    host  =  "https://wapps.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/"
    files  =  c(
        i450k_ns.xlsx  =  "48639-non-specific-probes-Illumina450k.xlsx",
        i450k_pl.xlsx  =  "48640-polymorphic-CpGs-Illumina450k.xlsx")
    
    for( i in seq_along(files)) {
        download.file(
            url  =  paste0(host, files[i]),
            destfile  =  paste0("DataDir/", names(files)[i]),
            mode  =  "wb",
            quiet  =  TRUE)
        }
    }

#. Ref Files ----
# library(readxl)
dataDir <- ("/Users/Emma/Documents/Bioinformatics/Epigenetics/DataDir/")
ex1  =  read_excel(paste(dataDir, "i450k_ns.xlsx", sep = "/"), sheet  =  1)
ex2  =  read_excel(paste(dataDir, "i450k_pl.xlsx", sep = "/"), sheet  =  1)
ex3  =  read_excel(paste(dataDir, "i450k_pl.xlsx", sep = "/"), sheet  =  2)

exclude.snp  =  unique(c(
    ex1$TargetID,
    ex2$PROBE,
    ex3$PROBE[ (ex3$BASE_FROM_SBE < 10) & (ex3$AF > 0.01)]))

# mSetSqFlt[[1]] == GSE60655, 450k
keep <- !(featureNames(mSetSqFlt[[1]]) %in% exclude.snp)
table(keep)
# FALSE   TRUE 
# 82024 385220
mSetSqFlt[[1]] <- mSetSqFlt[[1]][keep,] 
dim(mSetSqFlt[[1]])  # 385220     14
rm(host, files, i, ex1, ex2, ex3)
