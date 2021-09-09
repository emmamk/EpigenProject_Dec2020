#### pData ####
setwd('/Users/Emma/Documents/Bioinformatics/EpigenProject_Dec2020/GSE60655_End')

# library(GEOquery)
# geoMat <- getGEO("GSE60655", getGPL = FALSE)
# pD.all <- pData(geoMat[[1]])
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

pD <- pD.all[, c("title", "gender.ch1", "source_name_ch1", "group.ch1", 
                 "subject.ch1", "geo_accession", "slide.ch1", "array.ch1", "characteristics_ch1.3")]
head(pD, 3)
names(pD)[c(1:3,5,7:9)] <- c("rep", "gender", "timepoint", "sample", "slide", "array.position", "batch")

pD$rep <- str_extract(pD$rep, pattern="rep.")
pD$timepoint <- str_extract(pD$timepoint, pattern="before|after")
pD$timepoint <- factor(pD$timepoint, levels = c("before", "after"))
pD$sample <- paste0("S", pD$sample)
pD$ID <- paste(pD$sample, pD$timepoint, sep=".")
pD$batch <- str_replace(pD$batch, "batch: ", "")
pD$barcode <- paste(pD$slide, pD$array.position, sep = "_")
head(pD, 3)
#             rep gender timepoint group.ch1 sample geo_accession      slide array.position batch        ID           barcode
# GSM1484233 rep1   male    before        T1     S1    GSM1484233 5684819044         R01C01     A S1.before 5684819044_R01C01
# GSM1484234 rep1   male     after        T2     S1    GSM1484234 5684819044         R02C01     A  S1.after 5684819044_R02C01
# GSM1484235 rep1   male    before        T1     S4    GSM1484235 5684819044         R03C01     A S4.before 5684819044_R03C01



#### exploring the pdata ####
dim(pD)  # 36 11
table(pD$rep)
# rep1 rep2 
# 34    2

table(pD$slide)  # batch
# 5684819044 6164621122 6164621126 
# 12         12         12



#### subsetting the data ####
#.removing replicates ----
keep <- grep("rep1", pD$rep)
pD <- pD[keep,]; dim(pD)  # 34  11
table(pD$rep)
# rep1 
# 34


#. removing female subjects ----
keep <- grep("^male$", pD$gender)
pD <- pD[keep,]; dim(pD)  # 14 10
table(pD$gender)
# male 
# 14 


pD
#             rep gender timepoint group.ch1 sample geo_accession      slide array.position batch         ID           barcode
# GSM1484233 rep1   male    before        T1     S1    GSM1484233 5684819044         R01C01     A  S1.before 5684819044_R01C01
# GSM1484234 rep1   male     after        T2     S1    GSM1484234 5684819044         R02C01     A   S1.after 5684819044_R02C01
# GSM1484235 rep1   male    before        T1     S4    GSM1484235 5684819044         R03C01     A  S4.before 5684819044_R03C01
# GSM1484236 rep1   male     after        T2     S4    GSM1484236 5684819044         R04C01     A   S4.after 5684819044_R04C01
# GSM1484237 rep1   male    before        T1     S5    GSM1484237 5684819044         R05C01     A  S5.before 5684819044_R05C01
# GSM1484238 rep1   male     after        T2     S5    GSM1484238 5684819044         R06C01     A   S5.after 5684819044_R06C01
# GSM1484239 rep1   male    before        T1    S12    GSM1484239 5684819044         R01C02     A S12.before 5684819044_R01C02
# GSM1484240 rep1   male     after        T2    S12    GSM1484240 5684819044         R02C02     A  S12.after 5684819044_R02C02
# GSM1484241 rep1   male    before        T1    S13    GSM1484241 5684819044         R03C02     A S13.before 5684819044_R03C02
# GSM1484242 rep1   male     after        T2    S13    GSM1484242 5684819044         R04C02     A  S13.after 5684819044_R04C02
# GSM1484243 rep1   male    before        T1    S26    GSM1484243 5684819044         R05C02     A S26.before 5684819044_R05C02
# GSM1484244 rep1   male     after        T2    S26    GSM1484244 5684819044         R06C02     A  S26.after 5684819044_R06C02
# GSM1484245 rep1   male    before        T1     S7    GSM1484245 6164621122         R01C01     B  S7.before 6164621122_R01C01
# GSM1484246 rep1   male     after        T2     S7    GSM1484246 6164621122         R02C01     B   S7.after 6164621122_R02C01

