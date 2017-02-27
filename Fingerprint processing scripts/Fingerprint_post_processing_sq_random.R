# Post processing fingerprint_random datasets
# Author: Gabriel Altschuler
# Timestamp: 20110612
# Status: Completed 8th August 2011
# location: hpc111
#
# Processing the fingerprinted data
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# SCRIPT IS IN 3 DISTINCT PARTS
# quit and restart R between each to conserve memory
# It is not properly returned to system otherwise (even with garbage collection)

########
# PART 1
########

# First require matrix of raw single chip enrichments
# Using SCE processed data, mean rank method

# set directory holding fingerprint collection
fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/random_sq"
# define file identifier
header<-"random_sq_"
setwd(fingerpath)

# retrieve list of files, some will not have been properly processed
# use 200 bytes used as cutoff for successful fingerprint calculation

files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
unprocessedFiles<-dir(fingerpath)[file.info(dir(fingerpath))$size < 200]
files<-files[grep(header, files)]
unprocessedFiles<-unprocessedFiles[grep(header, unprocessedFiles)]
print(paste(length(files), "successfully fingerprinted GEO arrays"))
print(paste(length(unprocessedFiles), "GEO arrays could not be fingerprinted"))


############
# Remove arrays with <90% probe representation
library(GEOmetadb)
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

# retrieve metadata, including number of rows in the GEO data table
gsm.nrow = dbGetQuery(con, "select gsm, data_row_count from gsm")
gsm.gpl<-geoConvert(gsub(header, "", files), c("gpl"))$gpl
dbDisconnect(con)

probe.data<-as.data.frame(matrix(nrow = length(files), ncol = 3))
colnames(probe.data)<-c("GSM","GPL","nrow")
probe.data$GSM<-gsub(header, "", files)
probe.data$GPL<-gsm.gpl$to_acc[match(probe.data$GSM, gsm.gpl$from_acc)]
probe.data$nrow<-gsm.nrow$data_row_count[match(probe.data$GSM, gsm.nrow$gsm)]

# for each platform build a table of the data table nrows and find the max (i.e. complete array)
platforms<-levels(as.factor(probe.data$GPL))
platforms.nrow<-lapply(platforms, function(x){table(probe.data$nrow[probe.data$GPL == x])})
names(platforms.nrow)<-platforms
maxProbes<-sapply(platforms.nrow, function(x){max(as.numeric(names(x)))})

# define a TRUE/FALSE column for whether an array has at least 90% of the probes accounted for
probe.data$probes<-NA
for (i in 1:length(platforms)){
  threshold<-0.9*maxProbes[i]
  probe.data$probes[
    probe.data$GPL == platforms[i]
    ] <- (probe.data$nrow[
              probe.data$GPL == platforms[i]
              ] >= threshold)
}

# in addition, remove the Agilent arrays
# must remove Agilent arrays as these are two channel and cannot be analyzed in the same way
agilents<-c("GPL1708", "GPL4133", "GPL4134", "GPL6466", "GPL6480", "GPL6848", "GPL7202", "GPL887",  "GPL891")
probe.data$probes[probe.data$GPL %in% agilents]<-FALSE

# now re-define the file list and continue
files<-paste(header, probe.data$GSM[probe.data$probes], sep = "")
setwd(fingerpath)
length(files)
#######################
print("Normalizing using the POE approach")
print("constructing dataframes")
# construct list of fingerprints - SCGs (single chip gene set enrichments) - and platforms
# This displays a running total of the arrays processed for tracking
SCE<-vector("list", length(files))
platform<-vector("list", length(files))
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (i in 1:length(files)){
  load(files[i])
  print(i)
  # setTxtProgressBar(pb, i)
  SCE[[i]]<-temp1[[1]]$SCG
  platform[[i]]<-temp1[[1]]$platform
  }

# load geneset from data repository to give pathway names
geneset<-"random.v0.3.Hs.gs"
load(paste("/home/galtschu2/fingerprint/data/", geneset, sep = "")) 


# create SCG dataframe
names(SCE)<-gsub(header, "", files)
SCE.frame<-t(as.data.frame(SCE))
colnames(SCE.frame)<-names(get(geneset))
rownames(SCE.frame)<-gsub(header, "", files)
save(SCE.frame, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, ".frame.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

# create platform dataframe with timestamp
names(platform)<-gsub(header, "", files)
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

for (i in 1:ncol(SCE.frame)){
assign("pathway.SCE", SCE.frame[,i])
save(pathway.SCE, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/", header, 
                  ".pathway.SCE.", i, ".RData", sep = ""))
}


# Now quit R and restart to save on memory usage
# END OF PART 1 - UP TO HERE
###################

########
# PART 2
########

# this across 14 nodes
library(metaArray)
library(foreach)
library(multicore)
library(doMC)
registerDoMC(cores = 14)

# find the right files to load
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")

dir()[grep("platform", dir())]

#e.g.
load("platform.frame.2011-08-05.RData")

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
header<-"random_sq_"
# load list of pathways
geneset<-"random.v0.3.Hs.gs"
load(paste("/home/galtschu2/fingerprint/data/", geneset, sep = ""))
geneset<-get(geneset)

# create directory to store POE files
system(paste("mkdir ", "/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/POE/",
                                    header, sep = ""))

print("processing POE for each platform")
pb <- txtProgressBar(min = 0, max = length(geneset), style = 3)
for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("calculating POE for ", array, sep = ""))
  # for every pathway perform POE
  foreach (i = 1:length(geneset)) %dopar% {
    # load the SCG values for each pathway
    load(paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/",
                header, ".pathway.SCE.", i, ".RData", sep = ""))
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
    save(pathway.POE, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/POE/",
                                    header, "/pathway.SCE.", i, ".", array, ".RData", sep = ""))
    setTxtProgressBar(pb, i)            
    }
  }


######################
# Quit R and restart to conserve memory

# find the right files to load
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")

dir()[grep("platform", dir())]
load("platform.frame.2011-08-05.RData")

##########
# Convert platform names
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
##########

header<-"random_sq_"
# load list of pathways
geneset<-"random.v0.3.Hs.gs"
load(paste("/home/galtschu2/fingerprint/data/", geneset, sep = ""))
geneset<-get(geneset)

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
    load(paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/POE/",
                                    header, "/pathway.SCE.", i, ".", array, ".RData", sep = ""))
    POE.matrix[i,names(pathway.POE)]<-pathway.POE
    setTxtProgressBar(pb, i)
    }
  }

# remove arrays that are NA for the first pathway - this is usually due to a mis-match between the platform species and the annotated species
POE.matrix<-POE.matrix[,!(is.na(POE.matrix[1,]))]

save(POE.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, ".POE.matrix.",
                              Sys.Date(), ".RData",
                              sep = ""
                              )
    )


POE.matrix.0.001<-(POE.matrix > 0.001) - (POE.matrix < -0.001)
# how many NAs? - These are probably due to the species or chips where a pathway is not represented properly
(sum(is.na(POE.matrix.0.001)))/(dim(POE.matrix.0.001)[1]*dim(POE.matrix.0.001)[2])
# how many zeros?
(sum((POE.matrix.0.001 == 0), na.rm = TRUE))/(dim(POE.matrix.0.001)[1]*dim(POE.matrix.0.001)[2])

save(POE.matrix.0.001, file = paste("/data/shared/Fingerprint/", header, ".POE.matrix.0.001.",
                              Sys.Date(), ".RData",
                              sep = ""
                              )
    )

#################
