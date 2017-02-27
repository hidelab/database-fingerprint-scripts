# Normalizing fingerprint using POE_part2
#
# Author: Gabriel Altschuler
# Timestamp: 20110405
#
# to be run on hpc111
#
# Processing the fingerprinted data
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# This script had to be run in two parts to restore memory properly to the system

# run to calculate POE for each pathway, run in parallel using 14 cores
library(metaArray)
library(foreach)
library(multicore)
library(doMC)
registerDoMC(cores = 14)

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/")

# load list of platforms
load(grep("platform", dir(), value = TRUE))
platforms<-as.character(unique(platform.frame$Platform))

# load list of pathways
load("/home/galtschu2/fingerprint/data/kegg_wiki_TR_static.Hs.gs")
geneset<-kegg_wiki_TR_static.Hs.gs

print("processing POE for each platform")
pb <- txtProgressBar(min = 0, max = length(geneset), style = 3)
for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("calculating POE for ", array, sep = ""))
  # for every pathway perform POE
  foreach (i = 1:length(geneset)) %dopar% {
    # load the SCG values for each pathway
    load(paste("pathway.SCG.", i, ".RData", sep = ""))
    # subset for the array to be assessed
    sample <- pathway.SCG[grep(array, platform.frame$Platform)]
    # setup POE vector
    pathway.POE<-sample
    pathway.POE[]<-NA
    # subset to avoid NAs
    sample.valid<-sample[!(is.na(sample))]
    fit<-fit.em(sample.valid,
            cl = rep(0,length(sample.valid))
            )
    # write to POE file
    pathway.POE[!(is.na(sample))]<-fit$expr
    save(pathway.POE, file = paste("POE/pathway.SCG", i, array, "RData", sep = "."))
    setTxtProgressBar(pb, i)            
    }
  }


