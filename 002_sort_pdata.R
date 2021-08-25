#### Sorting pdata and rgSet ####
source("000_util.R")
source("001_data_download.R")
# library(stringr)
# library(dplyr)


#### GSE60655 ####
# pD.all
# library(minfi)
pD.GSE60655 <- as(pD.all[["GSE60655"]][, c("title", "gender:ch1", "source_name_ch1", "subject:ch1", "geo_accession")], "DataFrame")
# (as( ,"DataFrame"))が無いと後の pData(rgSet) <- pD がerror.
# pData(rgSet)はDataFrameである必要がある、とdeveloperが修正している:
# https://github.com/genomicsclass/labs/pull/90
# 参考: https://github.com/hansenlab/minfi/issues/174
head(pD.GSE60655)

# name: don't start with numbers nor contain "-"
names(pD.GSE60655)[c(1,2,3,4)] <- c("rep", "gender", "timepoint", "sample")
pD.GSE60655$rep <- str_extract(pD.GSE60655$rep, pattern="rep.")
pD.GSE60655$timepoint <- str_extract(pD.GSE60655$timepoint, pattern="before|after")
pD.GSE60655$sample <- paste0("S", pD.GSE60655$sample)
pD.GSE60655$ID <- paste(pD.GSE60655$sample, pD.GSE60655$timepoint, sep=".")

# removing female samples ----
pD.GSE60655 <- pD.GSE60655 %>% 
    subset(gender == "male")
with(pD.GSE60655, table(gender, timepoint))
#       after before
# male     8      8


head(pD.GSE60655, 3)
# DataFrame with 3 rows and 6 columns
#                    rep      gender   timepoint      sample geo_accession          ID
#            <character> <character> <character> <character>   <character> <character>
# GSM1484233        rep1        male      before          S1    GSM1484233   S1.before
# GSM1484234        rep1        male       after          S1    GSM1484234    S1.after
# GSM1484235        rep1        male      before          S4    GSM1484235   S4.before

tail(pD.GSE60655, 3)
# ataFrame with 3 rows and 6 columns
#                    rep      gender   timepoint      sample geo_accession          ID
#            <character> <character> <character> <character>   <character> <character>
# GSM1484246        rep1        male       after          S7    GSM1484246    S7.after
# GSM1484251        rep2        male      before         S13    GSM1484251  S13.before
# GSM1484252        rep2        male       after         S13    GSM1484252   S13.after



#### GSE114763 ####
pD.GSE114763 <- as(pD.all[["GSE114763"]][, c("title", "title", "geo_accession")], "DataFrame")
head(pD.GSE114763, 3)

# name: don't start with numbers nor contain "-"
names(pD.GSE114763)[c(1,2)] <- c("rep", "timepoint_org")
pD.GSE114763$rep <- tolower(str_extract(pD.GSE114763$rep, pattern="._Rep."))
pD.GSE114763$timepoint_org <- tolower(gsub("SkM_Epi_Mem_.*: Muscle_|7wk_|_P.*", "", pD.GSE114763$timepoint_org))
pD.GSE114763$sample <- paste0("S", str_sub(pD.GSE114763$rep, 1, 1))
pD.GSE114763$rep <- str_extract(pD.GSE114763$rep, pattern="rep.")

# selecting baseline and unloading for timepoint ----
pD.GSE114763 <- pD.GSE114763 %>% 
    subset(timepoint_org %in% c("baseline", "unloading"))
with(pD.GSE114763, table(sample, timepoint_org))
# timepoint
# sample baseline unloading
# S1        2         1      # S1 duplicated
# S2        1         1
# S3        1         1
# S4        1         1
# S5        1         1
# S6        1         1
# S7        1         1
# S8        1         1

# renaming timepoint
pD.GSE114763$timepoint <- ifelse(str_detect(pD.GSE114763$timepoint_org, "baseline"), "before", "after")
pD.GSE114763$ID <- paste(pD.GSE114763$sample, pD.GSE114763$timepoint, sep=".")

table(pD.GSE114763$rep)
# rep1 rep2 
# 16    1 

table(pD.GSE114763$sample)
# S1 S2 S3 S4 S5 S6 S7 S8 
# 3  2  2  2  2  2  2  2

table(pD.GSE114763$timepoint)
# after before
#     8      9

head(pD.GSE114763, 3)
# DataFrame with 3 rows and 5 columns
#                   rep   timepoint geo_accession      sample           ID
#            <character> <character>   <character> <character>  <character>
# GSM3149860        rep1      baseline    GSM3149860          S1      before   S1.before
# GSM3149863        rep1     unloading    GSM3149863          S1       after    S1.after
# GSM3149865        rep1      baseline    GSM3149865          S2      before   S2.before

tail(pD.GSE114763, 3)
# DataFrame with 3 rows and 5 columns
#                    rep   timepoint geo_accession      sample           ID
#            <character> <character>   <character> <character>  <character>
# GSM3149895        rep1      baseline    GSM3149895          S8      before   S8.before
# GSM3149898        rep1     unloading    GSM3149898          S8       after    S8.after
# GSM3149899        rep2      baseline    GSM3149899          S1      before   S1.before



#### Creating a list of pData ####
pD <- list(pD.GSE60655, pD.GSE114763)
names(pD) <- c("GSE60655", "GSE114763")
pD

# $GSE60655
# DataFrame with 16 rows and 6 columns
#                    rep      gender   timepoint      sample geo_accession          ID
#            <character> <character> <character> <character>   <character> <character>
# GSM1484233        rep1        male      before          S1    GSM1484233   S1.before
# GSM1484234        rep1        male       after          S1    GSM1484234    S1.after
# GSM1484235        rep1        male      before          S4    GSM1484235   S4.before
# GSM1484236        rep1        male       after          S4    GSM1484236    S4.after
# GSM1484237        rep1        male      before          S5    GSM1484237   S5.before
# ...                ...         ...         ...         ...           ...         ...
# GSM1484244        rep1        male       after         S26    GSM1484244   S26.after
# GSM1484245        rep1        male      before          S7    GSM1484245   S7.before
# GSM1484246        rep1        male       after          S7    GSM1484246    S7.after
# GSM1484251        rep2        male      before         S13    GSM1484251  S13.before
# GSM1484252        rep2        male       after         S13    GSM1484252   S13.after
# 
# $GSE114763
# DataFrame with 17 rows and 5 columns
#                    rep   timepoint geo_accession      sample           ID
#            <character> <character>   <character> <character>  <character>
# GSM3149860        rep1      baseline    GSM3149860          S1      before   S1.before
# GSM3149863        rep1     unloading    GSM3149863          S1       after    S1.after
# GSM3149865        rep1      baseline    GSM3149865          S2      before   S2.before
# GSM3149868        rep1     unloading    GSM3149868          S2       after    S2.after
# GSM3149870        rep1      baseline    GSM3149870          S3      before   S3.before
# ...                ...           ...           ...         ...         ...         ...
# GSM3149890        rep1      baseline    GSM3149890          S7      before   S7.before
# GSM3149893        rep1     unloading    GSM3149893          S7       after    S7.after
# GSM3149895        rep1      baseline    GSM3149895          S8      before   S8.before
# GSM3149898        rep1     unloading    GSM3149898          S8       after    S8.after
# GSM3149899        rep2      baseline    GSM3149899          S1      before   S1.before

names(pD)  # "GSE60655"  "GSE114763"



#### removing duplicated replication ####
for (i in seq_along(pD)) {
    pdata <- pD[[i]]
    rgset <- rgSet[[i]]
    cat(names(pD)[i], "\n")
    keep <- grep("rep1", pdata$rep)
    pdata <- pdata[keep,]
    rgset <- rgset[,keep]
    cat("table(pdata$rep):", table(pdata$rep), "\n")
    cat("dim(pD):", dim(pdata), "\n")  # 34  7
    cat("dim(rgSet):", dim(rgset), "\n\n")  # 622399     34
    
    sampleNames(rgset) <- rownames(pdata) <- pdata$ID
    pdata <- pdata[sampleNames(rgset),]
    pData(rgset) <- pdata

    pD[[i]] <- pdata
    rgSet[[i]] <- rgset
    
    cat("pData:", "\n")
    print(pD[[i]])
    cat("\n\n")
    
    cat("rgSet:", "\n")
    print(rgSet[[i]])
    cat("\n\n")
}

# memo:
# GSE60655 
# table(pdata$rep): 14 
# dim(pD): 14 6 
# dim(rgSet): 622399 14 
# 
# pData: 
# DataFrame with 14 rows and 6 columns
#                    rep      gender   timepoint      sample geo_accession          ID
#            <character> <character> <character> <character>   <character> <character>
# S1.before         rep1        male      before          S1    GSM1484233   S1.before
# S1.after          rep1        male       after          S1    GSM1484234    S1.after
# S4.before         rep1        male      before          S4    GSM1484235   S4.before
# S4.after          rep1        male       after          S4    GSM1484236    S4.after
# S5.before         rep1        male      before          S5    GSM1484237   S5.before
# ...                ...         ...         ...         ...           ...         ...
# S13.after         rep1        male       after         S13    GSM1484242   S13.after
# S26.before        rep1        male      before         S26    GSM1484243  S26.before
# S26.after         rep1        male       after         S26    GSM1484244   S26.after
# S7.before         rep1        male      before          S7    GSM1484245   S7.before
# S7.after          rep1        male       after          S7    GSM1484246    S7.after
# 
# rgSet: 
# class: RGChannelSet 
# dim: 622399 14 
# metadata(0):
# assays(2): Green Red
# rownames(622399): 10600313 10600322 ... 74810490 74810492
# rowData names(0):
# colnames(14): S1.before S1.after ... S7.before S7.after
# colData names(6): rep gender ... geo_accession ID
# Annotation
# array: IlluminaHumanMethylation450k
# annotation: ilmn12.hg19
# 
# 
# GSE114763 
# table(pdata$rep): 16 
# dim(pD): 16 5 
# dim(rgSet): 1051943 16
# 
# pData: 
# DataFrame with 16 rows and 6 columns
#                   rep timepoint_org geo_accession      sample   timepoint          ID
#           <character>   <character>   <character> <character> <character> <character>
# S1.before        rep1      baseline    GSM3149860          S1      before   S1.before
# S1.after         rep1     unloading    GSM3149863          S1       after    S1.after
# S2.before        rep1      baseline    GSM3149865          S2      before   S2.before
# S2.after         rep1     unloading    GSM3149868          S2       after    S2.after
# S3.before        rep1      baseline    GSM3149870          S3      before   S3.before
# ...               ...           ...           ...         ...         ...         ...
# S6.after         rep1     unloading    GSM3149888          S6       after    S6.after
# S7.before        rep1      baseline    GSM3149890          S7      before   S7.before
# S7.after         rep1     unloading    GSM3149893          S7       after    S7.after
# S8.before        rep1      baseline    GSM3149895          S8      before   S8.before
# S8.after         rep1     unloading    GSM3149898          S8       after    S8.after
# 
# rgSet:
# lass: RGChannelSet 
# dim: 1051943 16 
# metadata(0):
# assays(2): Green Red
# rownames(1051943): 1600101 1600111 ... 99810990 99810992
# rowData names(0):
# colnames(16): S1.before S1.after ... S8.before S8.after
# colData names(5): rep timepoint_org ... timepoint ID
# Annotation
# array: IlluminaHumanMethylationEPIC
# annotation: ilm10b4.hg19

