#### Preparing the data: GSE60655 & GSE114763 ####
setwd("/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/")
# source("000.util.R")


##### Downloading idat files from GEO #####
# library(GEOquery)
idatFiles <- list()
for (i in seq_along(gseid)) {
    cat("####", gsedir[i], "####", "\n")
    setwd(gsedir.path[i])
    gse <- gseid[i]
    
    if (dir.exists("idat")) {
        cat("idat files are already downloaded.", "\n")
        idatFiles[[gse]] <- list.files("idat", pattern = "idat$")
        cat("number of idat files:", length(idatFiles[[i]]), "\n\n") 
    } else {
        cat("downloading idat files....", "\n")
        getGEOSuppFiles(gse, makeDirectory = F)
        untar(paste0(gse, "_RAW.tar"), exdir = "idat")
        
        cat("finished downloading...", "\n")
        cat("decompressing idat.gz files...", "\n")
        idat.files <- list.files("idat", pattern = "idat.gz$", full = TRUE)
        sapply(idat.files, gunzip, overwrite = TRUE)
        
        idat.files <- list.files("idat", pattern = "idat")
        idatFiles[[gse]] <- idat.files
        cat("finished decompressing...", "\n")
        print(head(idat.files))
        cat("number of idat files:", length(idat.files), "\n\n")
        
    }
    
    if (i != length(gsedir)) { next }
    setwd(project.dir)
    
    }


#### Preparing idat files to use methylumIDAT (GSE60655) ####
if (sum(grepl("^GSM", idatFiles[["GSE60655"]])) > 0) {
    files <- list.files("GSE60655_End/idat", pattern = "idat", full = TRUE)
    files.renamed <- gsub("GSM[0-9]+_?", "", files)
    file.rename(files, files.renamed)
    idatFiles[["GSE60655"]] <- list.files("GSE60655_End/idat", pattern = "idat")
    rm(filse, files.renamed)
    }

barcodes <- str_extract(idatFiles[["GSE60655"]], "\\d+_R\\d+C\\d+"); length(barcodes) # 72    



#### Reading in idat files ####
mlset <- list()
for (i in seq_along(gseid)) {
    cat("#### reading in idat files:", gsedir[i], "####", "\n")
    setwd(gsedir.path[i])
    gse <- gseid[i]
    
    if (gse == "GSE60655") {
        mlset[[gse]] <- methylumIDAT(barcodes = barcodes, idatPath = "idat")
        cat("dim(mlset):", dim(mlset[[i]]), "\n")  # MethyLumiSet, 485577  36
        
        if (!file.exists("org.probes.gse60655.txt")) {
        write.table(featureNames(mlset), "./org.probes.gse60655.txt", quote = F, row.names = F, col.names = F)
        }
        
    } else { # GSE114763
        mlset[[gse]] <- readEPIC("idat", parallel = T, n = F)
        cat("dim(mlset) EPIC:", dim(mlset[[gse]]), "\n")  # MethyLumiSet, 866297  40
        
        #### Keeping the same probes with GSE60655 ####
        # the probes are the original
        keep.probes <- read.table("../GSE60655_End/org.probes.gse60655.txt")
        keep.probes <- keep.probes$V1
        keep <- featureNames(mlset[[gse]]) %in% keep.probes
        mlset[[gse]] <- mlset[[gse]][keep,]
        cat("dim(mlset) after adjusting to 450K:", dim(mlset[[gse]]), "\n\n")
        
    }
    
    cat("beta values:", "\n")
    print(mlset[[i]]@assayData$betas[1:3,1:3]); cat("\n")
    
    cat("detP values:", "\n")
    print(mlset[[i]]@assayData$pvals[1:3,1:3]); cat("\n\n")

    
    if (i != length(gsedir)) { next }
    setwd(project.dir)
    cat("finished creating mlsets.", "\n\n")
    
    }


mlset.raw <- mlset
mlset <- mlset.raw

