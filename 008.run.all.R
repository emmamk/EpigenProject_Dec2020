#### scripts for the project ####
project.dir <- "/Users/Emma/GitHub/EpigenProject_Dec2020/"
setwd(project.dir)
source("000.util.R")
source("001.data.prep.R")
source("002.pdata.sort.R")
source("003.quality.control.R")
source("004.normalizations.R")
source("005.differential.methylated.analysis.R")
source("006.gometh.go.overlaps.exclusives.R")
source("007.gometh.kegg.sig.overlaps.exclusives.R")

