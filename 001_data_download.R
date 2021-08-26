#### Scripts for Downloading Dataset ####
# workflow ref: https://kasperdanielhansen.github.io/genbioconductor/html/minfi.html
# source("000_util.R")


##### Obtaining the Data #####
# library(GEOquery)
# getGEOSuppFiles("GSE114763", makeDirectory = F) # EPIC
# untar("GSE114763_RAW.tar", exdir = "idat")
for (i in gsedir) {
    cat(i, "\n")
    print(head(list.files(paste0(i, "/idat"), pattern = "idat")))
    cat("number of files:", length(list.files(paste0(i, "/idat"), pattern = "idat")), "\n\n")
    }

idatFiles <- list()
for (i in gsedir) {
    cat(i, "\n")
    idatdir <- paste0(i, "/idat")
    idatFiles[[i]] <- list.files(idatdir, pattern = ".idat", full = TRUE)
    # idatFiles[[i]] <- list.files(idatdir, pattern = "idat.gz$", full = TRUE)
    # sapply(idatFiles[[i]], gunzip, overwrite = TRUE)
    print(head(idatFiles[[i]]))
    cat("number of files:", length(idatFiles[[i]]), "\n\n")
    }


##### Read Raw Data From IDAT Files ####
# (Ignore Warnings) 
# library(minfi)
rgSet <- list()
for (i in gsedir) {
    cat(i, "\n")
    idatdir <- paste0(i, "/idat")
    rgSet[[i]] <- read.metharray.exp(idatdir)
    print(rgSet[[i]])
    cat("\n")
    
    cat("sampleNames:", "\n")
    print(head(sampleNames(rgSet[[i]]), 2))
    cat("\n\n")
    }



##### Rename Samples with Descriptive Ones #####
# library(GEOquery)
geoMat <- list()
pD.all <- list()
if (getOption('timeout') == 60) {
    options(timeout = 500)
    }
getOption('timeout')  # [1] 500

for (i in gseid) {
    cat(i, "\n")
    geomat <- getGEO(i, getGPL = FALSE)
    geoMat[[i]] <- geomat
    pD.all[[i]] <- pData(geomat[[1]])
    cat("pData:", "\n")
    print(pD.all[[i]][1:5, 1:5])
    cat("\n\n")
    }

length(geoMat)  # 2
names(geoMat)  # "GSE60655"  "GSE114763"
length(pD.all)  # 2
names(pD.all)  # "GSE60655"  "GSE114763"

