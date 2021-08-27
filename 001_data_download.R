#### Scripts for Downloading Dataset ####
# workflow ref: https://kasperdanielhansen.github.io/genbioconductor/html/minfi.html
# source("000_util.R")


##### Obtaining the Data #####
# library(GEOquery)
# getGEOSuppFiles("GSE114763", makeDirectory = F) # EPIC
# untar("GSE114763_RAW.tar", exdir = "idat")
idatFiles <- list()
for (i in seq_along(gseid)) {
    cat("####", gsedir[i], "####", "\n")
    idatdir <- paste0(gsedir[i], "/idat")
    gse <- gseid[i]
    idatFiles[[gse]] <- list.files(idatdir, pattern = ".idat$", full = TRUE)
    # idatFiles[[i]] <- list.files(idatdir, pattern = "idat.gz$", full = TRUE)
    # sapply(idatFiles[[i]], gunzip, overwrite = TRUE)
    print(head(idatFiles[[i]]))
    cat("number of files:", length(idatFiles[[i]]), "\n\n")
    
    if (i != length(gsedir)) { next }
    cat("names(idatFiles):", names(idatFiles), "\n\n\n") 
    }



##### Read Raw Data From IDAT Files ####
# (Ignore Warnings) 
# library(minfi)
rgSet <- list()
for (i in seq_along(gseid)) {
    cat("####", gsedir[i], "####", "\n")
    idatdir <- paste0(gsedir[i], "/idat")
    gse <- gseid[i]
    rgSet[[gse]] <- read.metharray.exp(idatdir)
    print(rgSet[i])
    
    cat("sampleNames:", "\n")
    print(head(sampleNames(rgSet[[i]]), 2))
    cat("\n\n")
    
    if (i != length(gsedir)) { next }
    cat("names(rgSet):", names(rgSet), "\n\n\n")  # "GSE60655"  "GSE114763"
    }

#### GSE60655_End #### 
# $GSE60655
# class: RGChannelSet 
# dim: 622399 36 
# metadata(0):
# assays(2): Green Red
# rownames(622399): 10600313 10600322 ... 74810490 74810492
# rowData names(0):
# colnames(36): GSM1484233_5684819044_R01C01 GSM1484234_5684819044_R02C01 ... GSM1484267_6164621126_R05C02
# GSM1484268_6164621126_R06C02
# colData names(0):
# Annotation
# array: IlluminaHumanMethylation450k
# annotation: ilmn12.hg19
# 
# sampleNames: 
# [1] "GSM1484233_5684819044_R01C01" "GSM1484234_5684819044_R02C01"


#### GSE114763_Res #### 
# $GSE114763
# class: RGChannelSet 
# dim: 1051943 40 
# metadata(0):
# assays(2): Green Red
# rownames(1051943): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
# colnames(40): GSM3149860_201496860025_R01C01 GSM3149861_201496860025_R02C01 ...
# GSM3149898_201496860128_R07C01 GSM3149899_201496860128_R08C01
# colData names(0):
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19
# 
# sampleNames: 
# [1] "GSM3149860_201496860025_R01C01" "GSM3149861_201496860025_R02C01"



##### Rename Samples with Descriptive Ones #####
# library(GEOquery)
geoMat <- list()
pD.all <- list()
if (getOption('timeout') == 60) {
    options(timeout = 500)
    }
getOption('timeout')  # [1] 500

for (i in seq_along(gseid)) {
    cat("####", gsedir[i], "####", "\n")
    geomat <- getGEO(gseid[i], getGPL = FALSE)
    gse <- gseid[i]
    geoMat[[gse]] <- geomat
    pD.all[[gse]] <- pData(geomat[[1]])
    
    cat("####", gsedir[i], "####", "\n")
    cat("pData:", "\n")
    print(pD.all[[gse]][1:3,1:5]); cat("\n")
    
    if (i != length(gseid)) { next }
    
    cat("length(geoMat):", length(geoMat), "\n")  # 2
    cat("names(geoMat):", names(geoMat), "\n")  # "GSE60655"  "GSE114763"
    cat("length(pD.all):", length(pD.all), "\n")  # 2
    cat("names(pD.all):", names(pD.all), "\n\n\n")  # "GSE60655"  "GSE114763"
    }


#### GSE60655_End ####
# pData: 
#                       title geo_accession                status submission_date last_update_date
# GSM1484233 muscle_T1_1_rep1    GSM1484233 Public on Nov 04 2014     Aug 22 2014      Nov 05 2014
# GSM1484234 muscle_T2_1_rep1    GSM1484234 Public on Nov 04 2014     Aug 22 2014      Nov 05 2014
# GSM1484235 muscle_T1_4_rep1    GSM1484235 Public on Nov 04 2014     Aug 22 2014      Nov 05 2014


#### GSE114763_Res ####
# pData: 
#                                                             title geo_accession                status submission_date last_update_date
# GSM3149860    SkM_Epi_Mem_1: Muscle_Baseline_Participant_001_Rep1    GSM3149860 Public on May 23 2018     May 22 2018      May 23 2018
# GSM3149861       SkM_Epi_Mem_2: Muscle_Acute_Participant_001_Rep1    GSM3149861 Public on May 23 2018     May 22 2018      May 23 2018
# GSM3149862 SkM_Epi_Mem_3: Muscle_7wk_Loading_Participant_001_Rep1    GSM3149862 Public on May 23 2018     May 22 2018      May 23 2018
# 
# length(geoMat): 2 
# names(geoMat): GSE60655 GSE114763 
# length(pD.all): 2 
# names(pD.all): GSE60655 GSE114763 


