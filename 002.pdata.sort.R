#### Preparing pData ####
project.dir <- "/Users/Emma/GitHub/EpigenProject_Dec2020/"
setwd(project.dir)
# source("000.util.R")
# source("001.data.prep.R")

# library(GEOquery)
pD.all <- list()
pD <- list()
for (i in seq_along(gseid)) {
    cat("#### reading in pheno data:", gsedir[i], "####", "\n")
    gse <- gseid[i]
    sampleinfo <- paste(gsedir[[i]], "sampleInfo.txt", sep = "/")
    
    if (!file.exists(sampleinfo)) {
        geoMat <- getGEO(gse, getGPL = FALSE)
        pdall <- pData(geoMat[[1]])
        pD.all[[gse]] <- pdall
        write.table(pdall, sampleinfo, sep = "\t")
    } else {
        pdall <- read.table(sampleinfo, sep = "\t")
        pD.all[[gse]] <- pdall
        }
    
    # cat("treatment_protocol:", "\n")
    # cat(head(pdall$treatment_protocol_ch1, 1), "\n\n")
    
    # cat("data_processing column:", "\n")
    # cat(head(pdall$data_processing, 1), "\n\n")
    
    
    if (gse == "GSE60655") {
        pd <- pdall[, c("title", "gender.ch1", "source_name_ch1", "group.ch1", 
                         "subject.ch1", "geo_accession", "slide.ch1", "array.ch1", "characteristics_ch1.3")]
        names(pd)[c(1:3,5,7:9)] <- c("rep", "gender", "timepoint", "sample", "slide", "array.position", "batch")
        
        pd$rep <- str_extract(pd$rep, pattern="rep.")
        pd$timepoint <- str_extract(pd$timepoint, pattern="before|after")
        pd$timepoint <- factor(pd$timepoint, levels = c("before", "after"))
        pd$sample <- paste0("S", pd$sample)
        pd$ID <- paste(pd$sample, pd$timepoint, sep=".")
        pd$batch <- str_replace(pd$batch, "batch: ", "")
        pd$barcode <- paste(pd$slide, pd$array.position, sep = "_")
        pd$slide <- as.character(pd$slide)
        
        #. removing female subjects ----
        keep <- grep("^male$", pd$gender)
        pd <- pd[keep,]; dim(pd)  # 14 10
        
    } else {
        
        pd <- pdall[, c("title", "title", "supplementary_file", "geo_accession")]
        names(pd)[c(1:3)] <- c("rep", "timepoint", "file")
        
        pd$rep <- tolower(str_extract(pd$rep, pattern="._Rep."))
        pd$timepoint <- tolower(gsub("SkM_Epi_Mem_.*: Muscle_|7wk_|_P.*", "", pd$timepoint))
        pd$sample <- paste0("S", str_sub(pd$rep, 1, 1))
        pd$ID <- paste(pd$sample, pd$timepoint, sep=".")
        pd$file <- str_extract(pd$file, "GSM[0-9]*_.*idat")
        pd$barcode <- gsub("_...\\.idat", "", pd$file)
        pd$slide <- as.character(str_extract(pd$file, "2014[0-9]+"))
        pd$rep <- str_extract(pd$rep, "rep.")
        
        #. selecting the timepoints of interest ----
        keep <- pd$timepoint %in% c("baseline", "loading")
        pd <- pd[keep,]
        
        pd$timepoint <- ifelse(str_detect(pd$timepoint, "baseline"), "before", "after")
        pd$timepoint <- factor(pd$timepoint, levels = c("before", "after"))
        
    }
    
    
    #.removing replicates ----
    keep <- grep("rep1", pd$rep)
    pd <- pd[keep,]
    
    #### exploring the pdata ####
    print(str(pd)); cat("\n\n")
    
    pD[[gse]] <- pd
    
}
        


#### Adjusting the order of rownames(pD) & featureNames(mlset) ####
for (i in seq_along(gseid)) {
    cat("#### merging mlset and pData:", gsedir[i], "####", "\n")
    gse <- gseid[i]
    pd <- pD[[i]]
    
    #. matching orders of mlset & pData ----
    mlset[[gse]] <- mlset[[gse]][,match(pd$barcode, sampleNames(mlset[[gse]]))]
    pData(mlset[[gse]]) <- pd
    
    #. renaming rownames(pd) & sampleNames(mlset[[gse]]) ----
    rownames(pd) <- pd$ID
    sampleNames(mlset[[gse]]) <- rownames(pd)
    pD[[gse]] <- pd
    
    print(mlset[[gse]]); cat("\n\n")

    }


cat("\n\n")