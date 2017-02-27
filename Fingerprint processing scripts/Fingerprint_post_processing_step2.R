# Post processing fingerprint_step 2
# Author: Gabriel Altschuler
# Timestamp: 20111130
# Status: Re-compiled to be sourced
# location: hpc111
#
# Processing the fingerprinted data
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# SCRIPT IS IN 3 DISTINCT PARTS, THIS IS 2 OF 3
# quit and restart R between each to conserve memory
# It is not properly returned to system otherwise (even with garbage collection)
########
# PART 2
########


# define paths and headers
header<-"sq_"
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")
# Define the location of the pathway frame
# If this script is run on same day as part 1  the line below loads the correct file
# If day is different, find the right files to load and adjust accordingly
# load(paste("platform.frame", Sys.Date(), "RData", sep = "."))
dir()[grep("platform", dir())]
filename <- readline(prompt = "Enter platform filename ") 
print("Loading platform data")
load(filename)

# create directory to store POE files
Pathwaypath <- "/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/"
POEpath <- "/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/POE/"
try(system(paste("mkdir ", POEpath, header, sep = "")))

# this script is run across 14 nodes

library(metaArray)
library(foreach)
library(multicore)
library(doMC)
registerDoMC(cores = 14)


# easiest way to combine equivalent chips is just to create a new chiptype, e.g. GPLcombination

print("combining equivalent platforms")
platform.frame$Platform<-gsub("GPL96", "HGU133A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL571", "HGU133A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL3921", "HGU133A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL4685", "HGU133A", platform.frame$Platform)

platform.frame$Platform<-gsub("GPL4133", "G4112", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL6480", "G4112", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL1708", "G4112", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL6848", "G4112", platform.frame$Platform)

platform.frame$Platform<-gsub("GPL891", "G4121A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL6466", "G4121A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL4134", "G4122F", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL7202", "G4122F", platform.frame$Platform)

print(table(platform.frame$Platform))
platforms<-as.character(unique(platform.frame$Platform))

#####

# load list of pathways
library(pathprint)
pathwayNames<-names(get(genesets[[1]]))

print("processing POE for each platform")
pb <- txtProgressBar(min = 0, max = length(pathwayNames), style = 3)
for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("calculating POE for ", array, sep = ""))
  # for every pathway perform POE
  foreach (i = 1:length(pathwayNames)) %dopar% {
    # load the SCG values for each pathway
    load(paste(Pathwaypath, header, ".pathway.SCE.", i, ".RData", sep = ""))
    # subset for the array to be assessed
    sample <- pathway.SCE[grep(array, platform.frame$Platform)]
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
    save(pathway.POE, file = paste(POEpath, header, "/pathway.SCE.", i, ".", array, ".RData", sep = ""))
    setTxtProgressBar(pb, i)            
    }
  }


######################
# Quit R and restart to conserve memory