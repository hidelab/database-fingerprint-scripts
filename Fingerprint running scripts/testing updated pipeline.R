# testing fingerprinting scripts
setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
# script to test that the updated pipeline gives the same results
# this is using the old chipframe, prior to using AILUN for all annotations
# Author: Gabriel Altschuler
# Timestamp: 20110609
# Status: Completed


# run on server
source("/home/galtschu2/fingerprint/scripts/gabriel functions.R")
setwd("/home/galtschu2/fingerprint/Fingerprint running scripts")
# need to load pathprint for annotation script but load updated chipframe on top
library(pathprint)
load("/home/galtschu2/fingerprint/data/chipframe.RData", .GlobalEnv)

load("/home/galtschu2/fingerprint/data/genesets.RData")
sapply(genesets, function(x){load(paste("/home/galtschu2/fingerprint/data/", x, sep = ""), .GlobalEnv)})

# test using a sample GSM on GPL570

temp<-geo2fingerprint(
              GSM = "GSM18424", 
              GEOthreshold = FALSE,
              GEOpath = "/home/galtschu2/Documents/Databases/GEOfiles/",
              geneset = "KEGG and Wikipathways and static",
              enrichmentMethod = "SCE",
              transformation = "rank",
              statistic = "mean",
              normalizedScore = FALSE,
              progressBar = FALSE
              )
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/SCE/SCE_KEGG_Wiki_static_GSM18424")
# check that they are the same
tail(temp[[1]])
tail(temp1$GSM18424[[1]])

# this is fine, N.B. obviously the KEGG pathways differ as they have been re-annotated.

