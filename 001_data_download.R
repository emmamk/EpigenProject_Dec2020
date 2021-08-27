#### Scripts for Downloading Dataset ####
# workflow ref: https://kasperdanielhansen.github.io/genbioconductor/html/minfi.html
# source("000_util.R")


##### Obtaining the Data #####
# library(GEOquery)
# getGEOSuppFiles("GSE114763", makeDirectory = F) # EPIC
# untar("GSE114763_RAW.tar", exdir = "idat")
for (i in gsedir) {
    cat(i, "\n")
    print(head(list.files(paste0(i, "/idat"), pattern = "idat")), 3)
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
names(idatFiles)



##### Read Raw Data From IDAT Files ####
# (Ignore Warnings) 
# library(minfi)
rgSet <- list()
for (i in gsedir) {
    cat("####", i, "####", "\n")
    idatdir <- paste0(i, "/idat")
    rgSet[[i]] <- read.metharray.exp(idatdir)
    print(rgSet[i])
    cat("\n")
    
    cat("sampleNames:", "\n")
    print(head(sampleNames(rgSet[[i]]), 2))
    cat("\n\n")
    
    if (i != gseid[length(gseid)]) { next }
    print(names(rgSet))  # "GSE60655_End"  "GSE114763_Res"
    }



##### Rename Samples with Descriptive Ones #####
# library(GEOquery)
geoMat <- list()
pD.all <- list()
if (getOption('timeout') == 60) {
    options(timeout = 500)
    }
getOption('timeout')  # [1] 500

for (i in seq_along(gsedir)) {
    cat("####", gsedir[i], "####", "\n")
    geomat <- getGEO(gseid[i], getGPL = FALSE)
    gse <- gseid[i]
    geoMat[[gse]] <- geomat
    pD.all[[gse]] <- pData(geomat[[1]])
    
    cat("####", gsedir[i], "####", "\n")
    cat("pData:", "\n")
    print(pD.all[[gse]][1:3,1:5])
    cat("\n\n")
    
    if (i != length(gseid)) { next }
    
    cat("length(geoMat):", length(geoMat), "\n")  # 2
    cat("names(geoMat):", names(geoMat), "\n")  # "GSE60655"  "GSE114763"
    cat("length(pD.all):", length(pD.all), "\n")  # 2
    cat("names(pD.all):", names(pD.all), "\n")  # "GSE60655"  "GSE114763"
    }


#### GSE60655_End ####
# pData: 
#                       title geo_accession                status submission_date last_update_date
# GSM1484233 muscle_T1_1_rep1    GSM1484233 Public on Nov 04 2014     Aug 22 2014      Nov 05 2014
# GSM1484234 muscle_T2_1_rep1    GSM1484234 Public on Nov 04 2014     Aug 22 2014      Nov 05 2014
# GSM1484235 muscle_T1_4_rep1    GSM1484235 Public on Nov 04 2014     Aug 22 2014      Nov 05 2014


#### GSE114763_Res ####
# pData: 
#                                                             title geo_accession                status submission_date
# GSM3149860    SkM_Epi_Mem_1: Muscle_Baseline_Participant_001_Rep1    GSM3149860 Public on May 23 2018     May 22 2018
# GSM3149861       SkM_Epi_Mem_2: Muscle_Acute_Participant_001_Rep1    GSM3149861 Public on May 23 2018     May 22 2018
# GSM3149862 SkM_Epi_Mem_3: Muscle_7wk_Loading_Participant_001_Rep1    GSM3149862 Public on May 23 2018     May 22 2018
# last_update_date
# GSM3149860      May 23 2018
# GSM3149861      May 23 2018
# GSM3149862      May 23 2018
# 
# 
# length(geoMat): 2 
# names(geoMat): GSE60655 GSE114763 
# length(pD.all): 2 
# names(pD.all): GSE60655 GSE114763 


