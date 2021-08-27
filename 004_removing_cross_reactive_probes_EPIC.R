#### Removing Cross Reactive Probes (EPIC) ####
#. Cross Reactive Probes ----
# (multiple places in the genome)
# Pidsley et al. 2016: EPIC
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1#Sec22

# source("000_util.R")
# source("003_quality_control.R")

#. Obtaining Reference Files for EPIC ----
# https://bioc.ism.ac.jp/packages/devel/bioc/vignettes/ramwas/inst/doc/RW5a_matrix.html#probes-with-snps-and-in-cross-reactive-regions
if (!file.exists("DataDir")) {
    dir.create("DataDir")
}
if (!file.exists("DataDir/S1_cross_reactive.csv")) {
    host = "https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/"
    files = c(
        S1_cross_reactive.csv     = '13059_2016_1066_MOESM1_ESM.csv',
        S4_snp_cpg.csv            = '13059_2016_1066_MOESM4_ESM.csv',
        S5_snp_base_extension.csv = '13059_2016_1066_MOESM5_ESM.csv',
        S6_snp_body.csv           = '13059_2016_1066_MOESM6_ESM.csv')
    
    for( i in seq_along(files) )
        download.file(
            url = paste0(host, files[i]),
            destfile = names(files)[i],
            mode = "wb",
            quiet = TRUE)
}

#. Ref Files ----
snpcpgs1 = read.csv(paste(dataDir,'S1_cross_reactive.csv', sep="/"), stringsAsFactors = FALSE)
snpcpgs4 = read.csv(paste(dataDir,'S4_snp_cpg.csv', sep="/"), stringsAsFactors = FALSE)
snpcpgs5 = read.csv(paste(dataDir,'S5_snp_base_extension.csv', sep="/"), stringsAsFactors = FALSE)
snpcpgs6 = read.csv(paste(dataDir,'S6_snp_body.csv', sep="/"), stringsAsFactors = FALSE)

xReactiveProbes = unique(c(
    snpcpgs1$X,
    snpcpgs4$PROBE,
    snpcpgs5$PROBE,
    snpcpgs6$PROBE[
        pmax(snpcpgs6$VARIANT_START - snpcpgs6$MAPINFO, 
             snpcpgs6$MAPINFO - snpcpgs6$VARIANT_END) < 10]))

# mSetSqFlt[[2]] == GSE114763, EPIC
id <- "GSE114763"
keep <- !(featureNames(mSetSqFlt[[id]]) %in% xReactiveProbes)
table(keep)
# TRUE 
# 776698
mSetSqFlt[[id]] <- mSetSqFlt[[id]][keep,] 
mSetSqFlt[[id]]

rm(snpcpgs1, snpcpgs4, snpcpgs5, snpcpgs6)
