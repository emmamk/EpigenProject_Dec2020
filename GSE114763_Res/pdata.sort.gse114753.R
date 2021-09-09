#### pData ####
setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/test/0902ref.before/GSE114763_Res')

# library(GEOquery)
# geoMat <- getGEO("GSE114763", getGPL = FALSE)
# pD.all <- pData(geoMat[[1]])
# write.table(pD.all, "sampleInfo.txt", sep = "\t")
pD.all <- read.table("sampleInfo.txt", sep = "\t")
head(pD.all)
tail(pD.all)

# head(pD.all$treatment_protocol_ch1, 1)
# Healthy, sedentary subjects performed supervised one-legged, knee-extension exercise training for three months.
# Skeletal muscle biopsies from the vastus lateralis were taken at rest, before and after the training period.

# head(pD.all$data_processing, 1)
# Raw intesities exported from GenomeStudio were imported into R.
# Probes on chromosomes X and Y were removed.
# Probes with a detection P > 0.01 exceeding 5% of the samples were also filtered out. 
# Color-bias adjustment and quantile normalization were performed as implemented in lumi.
# BMIQ correction was applied to correct for probe type bias.
# Finally, M-values were batch-corrected usign the ComBat function.

pD <- pD.all[, c("title", "title", "supplementary_file", "geo_accession")]
head(pD, 3)
names(pD)[c(1:3)] <- c("rep", "timepoint", "file")
pD$rep <- tolower(str_extract(pD$rep, pattern="._Rep."))
pD$timepoint <- tolower(gsub("SkM_Epi_Mem_.*: Muscle_|7wk_|_P.*", "", pD$timepoint))
pD$sample <- paste0("S", str_sub(pD$rep, 1, 1))
pD$ID <- paste(pD$sample, pD$timepoint, sep=".")
str_view(pD$file, "2014[0-9]+.*[0-9]")
pD$file <- str_extract(pD$file, "GSM[0-9]*_.*idat")
pD$file.id <- str_replace(pD$file, "_Grn.idat", "")
pD$barcode <- str_extract(pD$file, "2014[0-9]+.*[0-9]")
pD$slide <- str_extract(pD$file, "2014[0-9]+")
pD$rep <- str_extract(pD$rep, "rep.")
head(pD, 3)
#             rep timepoint                                    file geo_accession sample          ID        slide             barcode
# GSM3149860 rep1  baseline GSM3149860_201496860025_R01C01_Grn.idat    GSM3149860     S1 S1.baseline 201496860025 201496860025_R01C01
# GSM3149862 rep1   loading GSM3149862_201496860025_R03C01_Grn.idat    GSM3149862     S1  S1.loading 201496860025 201496860025_R03C01
# GSM3149865 rep1  baseline GSM3149865_201496860025_R06C01_Grn.idat    GSM3149865     S2 S2.baseline 201496860025 201496860025_R06C01

tail(pD, 3)
#             rep timepoint                                    file geo_accession sample          ID        slide             barcode
# GSM3149892 rep1   loading GSM3149892_201496860128_R01C01_Grn.idat    GSM3149892     S7  S7.loading 201496860128 201496860128_R01C01
# GSM3149895 rep1  baseline GSM3149895_201496860128_R04C01_Grn.idat    GSM3149895     S8 S8.baseline 201496860128 201496860128_R04C01
# GSM3149897 rep1   loading GSM3149897_201496860128_R06C01_Grn.idat    GSM3149897     S8  S8.loading 201496860128 201496860128_R06C01



#### exploring the pdata ####
dim(pD)  # 40  8
table(pD$rep)
# rep1 rep2 
#  39    1

table(pD$slide)  # batch
# 201465940053 201496850072 201496860025 201496860106 201496860128 
# 8            8            8            8            8



#### subsetting the data ####
#. removing replicates ----
keep <- grep("rep1", pD$rep)
pD <- pD[keep,]; dim(pD)  # 39  7
table(pD$rep)
# rep1 
#   39


#. selecting the timepoints of interest ----
table(pD$timepoint %in% c("baseline", "loading"))
keep <- pD$timepoint %in% c("baseline", "loading")
pD <- pD[keep,]
table(pD$timepoint)
# baseline  loading 
# 8        8


pD
#             rep timepoint                                    file geo_accession sample          ID        slide             barcode
# GSM3149860 rep1  baseline GSM3149860_201496860025_R01C01_Grn.idat    GSM3149860     S1 S1.baseline 201496860025 201496860025_R01C01
# GSM3149862 rep1   loading GSM3149862_201496860025_R03C01_Grn.idat    GSM3149862     S1  S1.loading 201496860025 201496860025_R03C01
# GSM3149865 rep1  baseline GSM3149865_201496860025_R06C01_Grn.idat    GSM3149865     S2 S2.baseline 201496860025 201496860025_R06C01
# GSM3149867 rep1   loading GSM3149867_201496860025_R08C01_Grn.idat    GSM3149867     S2  S2.loading 201496860025 201496860025_R08C01
# GSM3149870 rep1  baseline GSM3149870_201465940053_R03C01_Grn.idat    GSM3149870     S3 S3.baseline 201465940053 201465940053_R03C01
# GSM3149872 rep1   loading GSM3149872_201465940053_R05C01_Grn.idat    GSM3149872     S3  S3.loading 201465940053 201465940053_R05C01
# GSM3149875 rep1  baseline GSM3149875_201465940053_R08C01_Grn.idat    GSM3149875     S4 S4.baseline 201465940053 201465940053_R08C01
# GSM3149877 rep1   loading GSM3149877_201496850072_R02C01_Grn.idat    GSM3149877     S4  S4.loading 201496850072 201496850072_R02C01
# GSM3149880 rep1  baseline GSM3149880_201496850072_R05C01_Grn.idat    GSM3149880     S5 S5.baseline 201496850072 201496850072_R05C01
# GSM3149882 rep1   loading GSM3149882_201496850072_R07C01_Grn.idat    GSM3149882     S5  S5.loading 201496850072 201496850072_R07C01
# GSM3149885 rep1  baseline GSM3149885_201496860106_R02C01_Grn.idat    GSM3149885     S6 S6.baseline 201496860106 201496860106_R02C01
# GSM3149887 rep1   loading GSM3149887_201496860106_R04C01_Grn.idat    GSM3149887     S6  S6.loading 201496860106 201496860106_R04C01
# GSM3149890 rep1  baseline GSM3149890_201496860106_R07C01_Grn.idat    GSM3149890     S7 S7.baseline 201496860106 201496860106_R07C01
# GSM3149892 rep1   loading GSM3149892_201496860128_R01C01_Grn.idat    GSM3149892     S7  S7.loading 201496860128 201496860128_R01C01
# GSM3149895 rep1  baseline GSM3149895_201496860128_R04C01_Grn.idat    GSM3149895     S8 S8.baseline 201496860128 201496860128_R04C01
# GSM3149897 rep1   loading GSM3149897_201496860128_R06C01_Grn.idat    GSM3149897     S8  S8.loading 201496860128 201496860128_R06C01

table(pD$slide)
# 201465940053 201496850072 201496860025 201496860106 201496860128 
# 3            3            4            3            3 

