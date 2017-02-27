# Normalizing fingerprint using POE
#
# Author: Gabriel Altschuler
# Timestamp: 20110401
#
# to be run on hpc111
#
# Processing the fingerprinted data
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# First require matrix of raw single chip enrichments
# Using SCE processed data, mean rank method

# set directory holding fingerprint collection
fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/SCE/"
setwd(fingerpath)

# retrieve list of files, some will not have been properly processed
# use 200 bytes used as cutoff for successful fingerprint calculation


files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
unprocessedFiles<-dir(fingerpath)[file.info(dir(fingerpath))$size < 200]
files<-files[grep("SCE_KEGG_Wiki_static_", files)]
unprocessedFiles<-unprocessedFiles[grep("SCE_KEGG_Wiki_static_", unprocessedFiles)]

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
if (!(exists("kegg_wiki_TR_static.Hs.gs"))){
load("/home/galtschu2/fingerprint/data/kegg_wiki_TR_static.Hs.gs", .GlobalEnv)
}

# create SCG dataframe
names(SCG)<-gsub("KEGG_Wiki_static_", "", files)
SCG.frame<-t(as.data.frame(SCG))
colnames(SCG.frame)<-names(kegg_wiki_TR_static.Hs.gs)
rownames(SCG.frame)<-gsub("KEGG_Wiki_static_", "", files)
save(SCG.frame, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/SCG.frame",
    Sys.Date(), "RData",
    sep = "."
    )
  )


# create platform dataframe with timestamp
names(platform)<-gsub("KEGG_Wiki_static_", "", files)
platform.frame<-data.frame(GEO = names(platform), Platform = unlist(platform))
save(platform.frame, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/platform.frame",
                                  Sys.Date(), "RData",
                                  sep = "."
                                  )
                              )

# each platform needs to be processed separately


##################################
# Further processing - this is optional
# it is possible that certain chips can be combined
# for example, HGU133A and HGU133A2 and the HT versions HTHGU133A and HTHG133A2
# these correspond to GPL96, GPL571, GPL3921, GPL4685
# very tempting to combine these chips as the same genes are represented
# this is also probably possible for the feature and probe versions of the Agilent arrays
# GPL4133  Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Feature Number version)
# GPL6480  Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)
# GPL1708  Agilent-012391 Whole Human Genome Oligo Microarray G4112A (Feature Number version)
# GPL6848  Agilent-012391 Whole Human Genome Oligo Microarray G4112A (Probe Name version)
# GPL891  Agilent-011978 Mouse Microarray G4121A (Feature Number version)
# GPL6466  Agilent-011978 Mouse Microarray G4121A (Probe Name version)
# GPL4134  Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (Feature Number version)
# GPL7202  Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (Probe Name version)
# this may also be possible for the illumina arrays
# The twelve-sample HumanHT-12 BeadChip is designed with the same panel of probes as the HumanWG-6 v3.0 BeadChip


# what is the distribution of the numbers?

table(platform.frame$Platform)
platforms.standard<-as.character(unique(platform.frame$Platform))

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
#################################

# Process data using the POE routine within the metaarray package
library(metaArray)

# problem that there are are small number (10) fingerprints that are all NAs
# remove these
platform.frame<-platform.frame[!(is.na(SCG.frame[,1])),]
SCG.frame<-SCG.frame[!(is.na(SCG.frame[,1])),]

# now have two issues
# need to run in parallel but dataframe is too large to process
# split into constituent pathway parts
for (i in 1:ncol(SCG.frame)){
  assign("pathway.SCG", SCG.frame[,i])
  save(pathway.SCG, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/",
                    "pathway.SCG.", i, ".RData", sep = ""))
  }

# save the platform file into the same location
save(platform.frame, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/platform.frame",
    Sys.Date(), "RData",
    sep = "."
    )
)

#####
# Now quit R and restart
# clearing workspace does not fully clean the memory

rm(list = ls(all = TRUE))

###################
# everything below here did not run properly

# create data matrices for POE and for POE.Pi
SCG.POE<-matrix(nrow = nrow(SCG.frame), ncol = ncol(SCG.frame))
colnames(SCG.POE)<-colnames(SCG.frame)
rownames(SCG.POE)<-rownames(SCG.frame)

# run routine - run in parallel, just use 2 cores due to memory restrictions
library(foreach)
library(multicore)
library(doMC)
registerDoMC(cores = 2)

print("processing POE for each platform")
pb <- txtProgressBar(min = 0, max = ncol(SCG.frame), style = 3)
for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("calculating POE for ", array, sep = ""))
  foreach (i = 1:ncol(SCG.frame)) %dopar% {
    sample <- SCG.frame[grep(array, platform.frame$Platform), i]
    sample.valid<-sample[!(is.na(sample))]
    fit<-fit.em(t(sample.valid),
            cl = rep(0,length(sample.valid))
            )
    SCG.POE[grep(array, platform.frame$Platform)[!(is.na(sample))], i]<-t(fit$expr)
    setTxtProgressBar(pb, i)        
    }
  }

save(SCG.POE, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/SCG.POE",
                          Sys.Date(), "RData",
                          sep = "."
                          )
                      )

SCG.POE1[grep(array, platform.frame$Platform)[!(is.na(sample))], i]<-fit$expr


test<-vector("numeric", 10)

foreach (i = 1:5) %dopar% {print(i); a<-i+10; print(a); test[i]<-a; print(test)}
foreach (i = 1:5) %do% {print(i); a<-i+10; print(a); test[i]<-a; print(test)}

for (i in 1:5) {print(i); a<-i+10; print(a); test[i]<-a; print(test)}

testfun<-function(x){print(x); a<-x+10; print(x); test[x]<-a; print(test)}
environment(testfun)<-globalenv()
mclapply(1:5, testfun)