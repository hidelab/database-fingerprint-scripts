# Normalizing fingerprint using POE - experimental datasets
#
# Author: Gabriel Altschuler
# Timestamp: 20110415
#
# to be run on hpc111
#
# Processing the fingerprinted data
# Experimental datasets
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# First require matrix of raw single chip enrichments
# Using SCE processed data, mean rank method

# set directory holding fingerprint collection
fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/SCE/experimental"
setwd(fingerpath)

# retrieve list of files, some will not have been properly processed
# use 200 bytes used as cutoff for successful fingerprint calculation

files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
unprocessedFiles<-dir(fingerpath)[file.info(dir(fingerpath))$size < 200]
files<-files[grep("SCE_Experimental_", files)]
unprocessedFiles<-unprocessedFiles[grep("SCE_Experimental_", unprocessedFiles)]

print(paste(length(files), "successfully fingerprinted GEO arrays"))
print(paste(length(unprocessedFiles), "GEO arrays could not be fingerprinted"))

print("Normalizing using the POE approach")
print("constructing dataframes")
# construct list of fingerprints - SCGs (single chip gene set enrichments) - and platforms
SCG<-vector("list", length(files))
platform<-vector("list", length(files))
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (i in 1:length(files)){
  load(files[i])
  print(i)
  # setTxtProgressBar(pb, i)
  SCG[[i]]<-temp1[[1]]$SCG
  platform[[i]]<-temp1[[1]]$platform
  }

# load geneset from data repository to give pathway names
if (!(exists("ExperimentalGeneSigs.Hs.gs"))){
load("/home/galtschu2/fingerprint/data/ExperimentalGeneSigs.Hs.gs", .GlobalEnv)
}


# create SCG dataframe
names(SCG)<-gsub("SCE_Experimental_", "", files)
SCG.frame<-t(as.data.frame(SCG))
colnames(SCG.frame)<-names(ExperimentalGeneSigs.Hs.gs)
rownames(SCG.frame)<-gsub("SCE_Experimental_", "", files)
save(SCG.frame, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/SCG.Experimental.frame",
    Sys.Date(), "RData",
    sep = "."
    )
  )


# create platform dataframe with timestamp
names(platform)<-gsub("SCE_Experimental_", "", files)
platform.frame<-data.frame(GEO = names(platform), Platform = unlist(platform))
save(platform.frame, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/platform.frame",
                                  Sys.Date(), "RData",
                                  sep = "."
                                  )
                              )

####
# each platform needs to be processed separately
# Certain platforms are equivalent and can be combined
# easiest way to combine is just to create a new chiptype, e.g. GPLcombination

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

table(platform.frame$Platform)
platforms<-as.character(unique(platform.frame$Platform))
#####

library(metaArray)

for (i in 1:ncol(SCG.frame)){
assign("pathway.SCG", SCG.frame[,i])
save(pathway.SCG, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/Experimental.",
                  "pathway.SCG.", i, ".RData", sep = ""))
}


###################
# Now quit R and restart to save on memory usage
# should really have run this with 15 nodes as there are 15 pathways :(
library(metaArray)
library(foreach)
library(multicore)
library(doMC)
registerDoMC(cores = 14)

# find the right files to load
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")

dir()[grep("platform", dir())]
load("platform.frame.2011-04-15.RData")

# load(grep("Experimental", dir(), value = TRUE))
# easiest way to combine is just to create a new chiptype, e.g. GPLcombination

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

table(platform.frame$Platform)
platforms<-as.character(unique(platform.frame$Platform))
#####

# load list of pathways
load("/home/galtschu2/fingerprint/data/ExperimentalGeneSigs.Hs.gs")
geneset<-ExperimentalGeneSigs.Hs.gs

print("processing POE for each platform")
pb <- txtProgressBar(min = 0, max = length(geneset), style = 3)
for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("calculating POE for ", array, sep = ""))
  # for every pathway perform POE
  foreach (i = 1:length(geneset)) %dopar% {
    # load the SCG values for each pathway
    load(paste("Experimental.pathway.SCG.", i, ".RData", sep = ""))
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
    save(pathway.POE, file = paste("POE/Experimental/pathway.SCG", i, array, "RData", sep = "."))
    setTxtProgressBar(pb, i)            
    }
  }

######################
# Quit R and restart to conserve memory


# find the right files to load
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")

dir()[grep("platform", dir())]
load("platform.frame.2011-04-15.RData")

# load(grep("Experimental", dir(), value = TRUE))
# easiest way to combine is just to create a new chiptype, e.g. GPLcombination

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

table(platform.frame$Platform)
platforms<-as.character(unique(platform.frame$Platform))

# load list of pathways
load("/home/galtschu2/fingerprint/data/ExperimentalGeneSigs.Hs.gs")
geneset<-ExperimentalGeneSigs.Hs.gs

# setup POE.matrix
POE.matrix<-matrix(nrow = length(geneset), ncol = length(platform.frame$Platform))
rownames(POE.matrix)<-names(geneset)
colnames(POE.matrix)<-platform.frame$GEO

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/")

pb <- txtProgressBar(min = 0, max = length(geneset), style = 3)
for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("Inserting POE values for ", array, sep = ""))
  for (i in 1:length(geneset)){
    load(paste("POE/Experimental/pathway.SCG", i, array, "RData", sep = "."))
    POE.matrix[i,names(pathway.POE)]<-pathway.POE
    setTxtProgressBar(pb, i)
    }
  }

save(POE.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/Experimental.POE.matrix",
                              Sys.Date(), "RData",
                              sep = "."
                              )
    )

#################
# save thresholded matrix and POE matrix into shared directory
# use thresholds 0.01 and 0.001
POE.matrix.0.1<-(POE.matrix > 0.1) - (POE.matrix < -0.1)
POE.matrix.0.01<-(POE.matrix > 0.01) - (POE.matrix < -0.01)
POE.matrix.0.001<-(POE.matrix > 0.001) - (POE.matrix < -0.001)
POE.matrix.0.0001<-(POE.matrix > 0.0001) - (POE.matrix < -0.0001)

# how many zeros?
(sum((POE.matrix.0.01 == 0), na.rm = TRUE))/(15 * 160837)
(sum((POE.matrix.0.001 == 0), na.rm = TRUE))/(15 * 160837)
(sum((POE.matrix.0.0001 == 0), na.rm = TRUE))/(15 * 160837)

save(POE.matrix, file = "/data/shared/Fingerprint/Experimental.POE.matrix.2011-04-15")
save(POE.matrix.0.01, file = "/data/shared/Fingerprint/Experimental.POE.matrix.0.01.2011-04-15")
save(POE.matrix.0.001, file = "/data/shared/Fingerprint/Experimental.POE.matrix.0.001.2011-04-15")
save(POE.matrix.0.0001, file = "/data/shared/Fingerprint/Experimental.POE.matrix.0.0001.2011-04-15")
