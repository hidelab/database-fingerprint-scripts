# Post processing fingerprint_step 1
# Author: Gabriel Altschuler
# Timestamp: 20121117
# Status: Re-compiled to be sourced
# location: hpc111
#
# Processing the fingerprinted data
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# SCRIPT IS IN 3 DISTINCT PARTS, THIS IS 1 OF 3
# quit and restart R between each to conserve memory
# It is not properly returned to system otherwise (even with garbage collection)

########
# PART 1
########

# First require matrix of raw single chip enrichments
# Using SCE processed data, mean rank method, DPD genesets

# set directory holding fingerprint collection
fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/DPD/"
# define file identifier
header<-"DPD_"
setwd(fingerpath)

# retrieve list of files, some will not have been properly processed
# use 200 bytes used as cutoff for successful fingerprint calculation

print("Loading files")

files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
unprocessedFiles<-dir(fingerpath)[file.info(dir(fingerpath))$size < 200]
files<-files[grep(header, files)]
unprocessedFiles<-unprocessedFiles[grep(header, unprocessedFiles)]
print(paste(length(files), "successfully fingerprinted GEO arrays"))
print(paste(length(unprocessedFiles), "GEO arrays could not be fingerprinted"))


############
# Remove arrays with <90% probe representation
print("Removing arrays with <90% probe representation")
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
load(files[1])
pathways <- rownames(temp1[[1]]$SCG)
print("Creating expression matrix")
# pre-allocate memory for matrix
SCE.matrix <- matrix(nrow = length(files), ncol = length(pathways))
colnames(SCE.matrix)<-pathways
rownames(SCE.matrix)<-gsub(header, "", files)
print(object.size(SCE.matrix), units = "auto")
SCE.matrix[,] <- 1
pb <- txtProgressBar(min = 0, max = length(files), style = 3)
for (i in 1:length(files)){
  load(files[i])
  print(i)
  # setTxtProgressBar(pb, i)
  SCE.matrix[i,]<-as.integer(temp1[[1]]$SCG)
  if (i %% 100 == 0){
  	gc()
  	}
  }
print("Saving matrix. The file size is")
print(object.size(SCE.matrix), units = "auto")
save(SCE.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, ".matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

# remove na rows
SCE.matrix <- SCE.matrix[,-c(1238, 18991)]
SCE.matrix <- SCE.matrix[-c(1718,6508),]
# calculate correlation matrix
# create line by line
# multiply by 10000000 to maintain integer values

cor.matrix <- matrix(nrow = ncol(SCE.matrix),
					 ncol = ncol(SCE.matrix))
colnames(cor.matrix)<-colnames(SCE.matrix)
rownames(cor.matrix)<-colnames(SCE.matrix)
print(object.size(cor.matrix), units = "auto")

imax = floor(ncol(SCE.matrix)/20) + 1		 
for (i in 1:imax){
  print(i)
  istart = 1 + 20 * (i-1)
  iend = istart + 19
  if (iend > ncol(SCE.matrix)){
  	iend = ncol(SCE.matrix)
  	}
  # setTxtProgressBar(pb, i)
  cor.matrix[,istart:iend]<-as.integer(1000000*cor(
  						SCE.matrix,
  						SCE.matrix[,istart:iend]),
  						use = "all.obs")
  	gc()
  }					 

save(cor.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, ".cor.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )


# perhaps for a start just try to produce a correlation matrix for the disease and hide genesets.
# can use this to look for the top anti-correlated drugs and signatures.

diseaseRef <- c(grep("disease", colnames(SCE.matrix)),
				  grep("Hide", colnames(SCE.matrix))
					)

# 700 pathways

disease.cor.matrix <- matrix(nrow = ncol(SCE.matrix),
					 ncol = 700)
colnames(disease.cor.matrix)<-colnames(SCE.matrix)[diseaseRef]
rownames(disease.cor.matrix)<-colnames(SCE.matrix)
print(object.size(disease.cor.matrix), units = "auto")

for (i in 1:35){
  print(i)
  istart = 1 + 20 * (i-1)
  iend = istart + 19
  refs = diseaseRef[istart:iend]
  # setTxtProgressBar(pb, i)
  system.time(
  disease.cor.matrix[,istart:iend]<-as.integer(1000000*cor(
  						SCE.matrix,
  						SCE.matrix[,refs]), use = "all.obs")
  						)
  	gc()
  }					 

save(disease.cor.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, "disease.cor.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

# also save a pathways cor matrix

pathwaysRef <- grep("Pathway\\.", colnames(cor.matrix))

# 633 pathways

pathway.cor.matrix <- cor.matrix[, pathwaysRef]
			 

save(pathway.cor.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, "pathway.cor.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )


# old script
#print("Normalizing using the POE approach")
#print("constructing dataframes")
## construct list of fingerprints - SCGs (single chip gene set enrichments) - and platforms
## This displays a running total of the arrays processed for tracking
#SCE<-vector("list", length(files))
#platform<-vector("list", length(files))
#pb <- txtProgressBar(min = 0, max = length(files), style = 3)
#for (i in 1:length(files)){
#  load(files[i])
#  print(i)
#  # setTxtProgressBar(pb, i)
#  SCE[[i]]<-temp1[[1]]$SCG
#  platform[[i]]<-temp1[[1]]$platform
#  }
#
## load geneset from data repository to give pathway names
#library(pathprint)
#
## create SCE dataframe
#print(paste("saving SCE dataframe with timestamp ", Sys.Date(), sep = " "))
#names(SCE)<-gsub(header, "", files)
#SCE.frame<-t(as.data.frame(SCE))
#colnames(SCE.frame)<-names(get(genesets[[1]]))
#rownames(SCE.frame)<-gsub(header, "", files)
#save(SCE.frame, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, ".frame.",
#    Sys.Date(), ".RData",
#    sep = ""
#    )
#  )

# remove as only one platform
## create platform dataframe with timestamp
#print(paste("saving platform dataframe with timestamp ", Sys.Date(), sep = " "))
#names(platform)<-gsub(header, "", files)
#platform.frame<-data.frame(GEO = names(platform), Platform = unlist(platform))
#save(platform.frame, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/platform.frame",
#                                  Sys.Date(), "RData",
#                                  sep = "."
#                                  )
#                              )

#####
#print("combining equivalent platforms")
## each platform needs to be processed separately
## Certain platforms are equivalent and can be combined
## easiest way to combine is just to create a new chiptype, e.g. GPLcombination
#
#platform.frame$Platform<-gsub("GPL96", "HGU133A", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL571", "HGU133A", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL3921", "HGU133A", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL4685", "HGU133A", platform.frame$Platform)
#
#platform.frame$Platform<-gsub("GPL4133", "G4112", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL6480", "G4112", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL1708", "G4112", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL6848", "G4112", platform.frame$Platform)
#
#platform.frame$Platform<-gsub("GPL891", "G4121A", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL6466", "G4121A", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL4134", "G4122F", platform.frame$Platform)
#platform.frame$Platform<-gsub("GPL7202", "G4122F", platform.frame$Platform)
#
#table(platform.frame$Platform)
#platforms<-as.character(unique(platform.frame$Platform))
######
#print("Saving pathway-specific data")
#library(metaArray)
#
#for (i in 1:ncol(SCE.frame)){
#assign("pathway.SCE", SCE.frame[,i])
#save(pathway.SCE, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/", header, 
#                  ".pathway.SCE.", i, ".RData", sep = ""))
#}
#

# Now quit R and restart to save on memory usage
# END OF PART 1
###################

