# Script for adding additional array types to the fingerprinting pipeline
# Author: Gabriel Altschuler
# Timestamp: 20111127
# Status: Updated for adding AFFY ST 1.0 platforms
# The file on the server that contains the chip information is 
# "/home/galtschu2/fingerprint/data/chipframe.RData"
# this is loaded into the global environment by the parent fingerprint script

# N.B. if the file name is changed this needs to be reflected in;
# GSEA working script within gabriel functions.R
# Reading in samples_server_update.R
# GEOsurvey_server_update.R

# The workflow has now been updated to solely use the AILUN annotation database

# load current chipframe
library(pathprint)
platforms<-names(chipframe)
# platforms to be added are the Human Gene 1.0 ST Array and Mouse Gene 1.0 ST Array  
platforms.new<-c("GPL6244", "GPL6246")
titles.new<-c("[HuGene-1_0-st] Affymetrix Human Gene 1.0 ST Array [transcript (gene) version]",
			  "[MoGene-1_0-st] Affymetrix Mouse Gene 1.0 ST Array [transcript (gene) version]")

# check they don't exist already
platforms.new %in% platforms

# The AILUN database has been downloaded into a local directory on hpc111
AILUN.dir<-"/data/shared/AILUN/"



chipframe.new<-vector("list", 0)
for (i in 1:length(platforms.new)){
   print(platforms.new[i])
   conn <- gzfile(description = paste(AILUN.dir, platforms.new[i], ".annot.gz", sep = ""), open = "rt")
   data <- read.delim(conn, stringsAsFactors = FALSE, header = FALSE)
   close(conn)
   annotation<-data.frame(EntrezID = data[,2], ID = data[,1], stringsAsFactors = FALSE)
   # update EntrezIDs - not totally necessary but highly advised
   # see entrezUpdate.R in GMAfunctions for details
   library(GMAfunctions)
   historyFile <- "hpc111"
   annotation$EntrezID<-entrezUpdate(annotation$EntrezID, historyFile = historyFile)
   # remove duplicated 
   annotation<-annotation[!duplicated(annotation$ID),]
   rownames(annotation)<-annotation$ID
   chipframe.new <- append(chipframe.new, list(list(ann = annotation, title = titles.new[i])))
}
names(chipframe.new)<-platforms.new



# add to existing chipframe
chipframe <- append(chipframe, chipframe.new)

# save into fingerprint repository directory
save(chipframe, file = "/home/galtschu2/fingerprint/data/chipframe.RData")

# this now needs to be built into a new version of the fingerprint R package to be used in the fingerprinting pipeline


# End of script
# previous instances below

#########################
# Previously used to update chipframe file from non-AILUN annotations 
##load("/home/galtschu2/fingerprint/data/chipframe.RData")
##save(chipframe, file = "/home/galtschu2/fingerprint/data/chipframe_pre_AILUN.RData")
source("/home/galtschu2/fingerprint/scripts/gabriel functions.R")
AILUN.dir<-"/home/galtschu2/Documents/Databases/AILUN/"
historyFile <- "hpc111"
platforms<-names(chipframe)
chipframe<-vector("list", 0)
for (i in 1:length(platforms)){
   print(platforms[i])
   conn <- gzfile(description = paste(AILUN.dir, platforms[i], ".annot.gz", sep = ""), open = "rt")
   data <- read.delim(conn, stringsAsFactors = FALSE, header = FALSE)
   close(conn)
   annotation<-data.frame(EntrezID = data[,2], ID = data[,1], stringsAsFactors = FALSE)
   annotation$EntrezID<-entrezUpdate(annotation$EntrezID, historyFile = historyFile)
   # remove duplicated 
   annotation<-annotation[!duplicated(annotation$ID),]
   rownames(annotation)<-annotation$ID
   chipframe <- append(chipframe, list(list(ann = annotation)))
}
names(chipframe)<-platforms

save(chipframe, file = "/home/galtschu2/fingerprint/data/chipframe.RData")


# end of script

#########################################################################################

###########
#OLD WORKFLOW
###########
# this is how the current dataframe was constructed
GPL570.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU133Plus2_Hs_ENTREZG_mapping.txt")
GPL1261.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/Mouse4302_Mm_ENTREZG_mapping.txt")
GPL339.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/MOE430A_Mm_ENTREZG_mapping.txt")
GPL96.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU133A_Hs_ENTREZG_mapping.txt")
GPL97.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU133B_Hs_ENTREZG_mapping.txt")
GPL81.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/MGU74Av2_Mm_ENTREZG_mapping.txt")
GPL8321.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/Mouse430A2_Mm_ENTREZG_mapping.txt")
GPL8300.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU95Av2_Hs_ENTREZG_mapping.txt")
GPL571.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU133A2_Hs_ENTREZG_mapping.txt")
GPL340.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/MOE430B_Mm_ENTREZG_mapping.txt")
load("/home/galtschu2/Documents/Databases/customCDFmappings/GPL2986.ann.RData")
chiplist<-c("GPL570", "GPL1261", "GPL339", "GPL96", "GPL97", "GPL81", "GPL8321", "GPL8300", "GPL571", "GPL340", "GPL2986")

annotation<-list(
GPL570.ann<-annotation(GPL570.data),
GPL1261.ann<-annotation(GPL1261.data),
GPL339.ann<-annotation(GPL339.data),
GPL96.ann<-annotation(GPL96.data),
GPL97.ann<-annotation(GPL97.data),
GPL81.ann<-annotation(GPL81.data),
GPL8321.ann<-annotation(GPL8321.data),
GPL8300.ann<-annotation(GPL8300.data),
GPL571.ann<-annotation(GPL571.data),
GPL340.ann<-annotation(GPL340.data),
GPL2986.ann<-GPL2986.ann
)
names(annotation) <- chiplist


chipframe = vector("list", length = 0)
for (i in 1:length(chiplist)){
  chipframe <- append(chipframe, list(list(ann = annotation[[chiplist[i]]])))
	}
names(chipframe)<-chiplist

save(chipframe, file = "data/chipframe_update.RData")




load("/Users/GabrielAltschuler/Documents/Databases/Gene sets/chipframe_update.RData")

# this is a list of dataframes, colnames(c("EntrezID", "AFFYID"))
# The AFFYID column is a list of all the probe ids with e.g. "1007_s_at","1053_at"
# EntrezID contains EntrezIDs, also suffixed with _at
# This is used to aggregate a GEO table using the script customCDFAnn
# The _at is stripped out before being passed on

# To insert new array need to consider
# 1) how the data is presented from GEOquery
# 2) how the platform is curated in GEO

# use custom CDF chips for Affy
# no custom CDF chips for Illumina arrays
# Use GEO annotations for Illumina arrays
# The required GEO platform table colnames are ID and ENTREZ_GENE_ID
# example below (N.B. takes a while to run)

# blank spaces indicate no match, need to remove these, added script platform2Ann to do this


# Check that this is the case for the Illumina arrays
# adding the following arrays
# GPL6947  Illumina HumanHT-12 V3.0 expression beadchip
# GPL6883  Illumina HumanRef-8 v3.0 expression beadchip
# GPL6104	Illumina humanRef-8 v2.0 expression beadchip
# GPL6102  Illumina human-6 v2.0 expression beadchip
# GPL6884  Illumina HumanWG-6 v3.0 expression beadchip

# GPL6887  Illumina MouseWG-6 v2.0 expression beadchip
# GPL6885  Illumina MouseRef-8 v2.0 expression beadchip
# GPL6103  Illumina mouseRef-8 v1.1 expression beadchip
# GPL6105  Illumina mouse-6 v1.1 expression beadchip


# test on server
source("/home/galtschu2/fingerprint/scripts/gabriel functions.R")
test<-platform2Ann("GPL571", onServer = TRUE)
test1<-platform2Ann("GPL6947", onServer = TRUE)
test2<-platform2Ann("GPL4133", onServer = TRUE)
test3<-platform2Ann("GPL6480", onServer = TRUE)

# script breaks down on second one, the Illumina array - must be lack of correct colnames
# changed so searches for columns called ID and *entrez* (not case sensitive)

# okay, so now have additional annotation files
# now need to add these to the chipframe
# rename chipframe_update_23_3_10.RData
# slight difference in that factors are not used for the EntrezID column
# Shouldn't make a difference, would prefer as characters so won't change for now

load("/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update.RData")
chiplist<-names(chipframe)

illuminaList<-c("GPL6947", "GPL6883", "GPL6104", "GPL6102", "GPL6884", "GPL6887", "GPL6885", "GPL6103", "GPL6105")

for (i in 1:length(illuminaList)){
  chipframe <- append(chipframe, list(list(ann = platform2Ann(illuminaList[i], onServer = TRUE))))
  }
names(chipframe)<-c(chiplist, illuminaList)

save(chipframe, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update_23_3_10.RData")

############

# Adding agilent arrays

# 29_3_11 adding additional Agilent arrays


# GPL4133  Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Feature Number version)
# GPL6480  Agilent-014850 Whole Human Genome Microarray 4x44K G4112F (Probe Name version)
# GPL1708	Agilent-012391 Whole Human Genome Oligo Microarray G4112A (Feature Number version)
# GPL6848  Agilent-012391 Whole Human Genome Oligo Microarray G4112A (Probe Name version)
# GPL887  Agilent-012097 Human 1A Microarray (V2) G4110B (Feature Number version)
# GPL891  Agilent-011978 Mouse Microarray G4121A (Feature Number version)
# GPL6466  Agilent-011978 Mouse Microarray G4121A (Probe Name version)
# GPL4134  Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (Feature Number version)
# GPL7202  Agilent-014868 Whole Mouse Genome Microarray 4x44K G4122F (Probe Name version)

# may be able to combine the Agilent probe and feature versions at the final processing stage
# test on server
source("/home/galtschu2/Documents/Databases/gabriel functions.R")
test<-platform2Ann("GPL571", onServer = TRUE)
test1<-platform2Ann("GPL6947", onServer = TRUE)
test2<-platform2Ann("GPL4133", onServer = TRUE)
test3<-platform2Ann("GPL6480", onServer = TRUE)

load("/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update_23_3_10.RData")

chiplist<-names(chipframe)
agilentList<-c("GPL4133", "GPL6480", "GPL1708", "GPL6848", "GPL887", "GPL891", "GPL6466", "GPL4134", "GPL7202")

for (i in 1:length(agilentList)){
  chipframe <- append(chipframe, list(list(ann = platform2Ann(agilentList[i], onServer = TRUE))))
  }
names(chipframe)<-c(chiplist, agilentList)

save(chipframe, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update_29_3_10.RData")

# want to add the high throughput chips, HTHGU133A and HTHGU133
# GPL3921 [HT_HG-U133A] Affymetrix HT Human Genome U133A Array
# GPL4685 HT-HG_U133A [U133AAofAv2] Affymetrix GeneChip HT-HG_U133A Early Access Array
# These are the same Affymetrix part number 
# noted that hgu133a, a2 and ht all have the same probe to entrez gene id matchings
# While the separate cdfs must be used to create the gene files, it should be fine just to use one for the probe to gene mapping

# also stop versioning these chipframes, just store one as chipframe.RData, with backup as timestamped
load(file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update_29_3_10.RData")

chipframe$GPL3921<-chipframe$GPL96
chipframe$GPL4685<-chipframe$GPL96

save(chipframe, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update_31_3_10.RData")
save(chipframe, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe.RData")



##########
# new session
# adding additional platforms 7th June 2011, including model organisms
library(pathprint)
data("chipframe")
# N.B. removed Sus scrofa from list as do not have homologene information
newPlatforms<-c("GPL1319", "GPL200", "GPL72", "GPL1322", "GPL341", "GPL85", "GPL1355", "GPL2700", "GPL2995", "GPL6333", "GPL91")
newPlatforms.list<-vector("list", 0)
for (i in 1:length(newPlatforms)){
  newPlatforms.list <- append(newPlatfoms.list, list(list(ann = platform2Ann(newPlatforms[i], onServer = TRUE))))
  }
names(newPlatforms.list)<-newPlatforms

##########
# New strategy - use the AILUN derived annotation files
# only problem is that the probe IDs are not unique
# Resolve by taking the first one

# ftp ftp://ailun.stanford.edu
# login with anonymous, anonymous
# cd ailun/annotation/geo/
# prompt n
# mget *
source("/home/galtschu2/fingerprint/scripts/gabriel functions.R")
AILUN.dir<-"/home/galtschu2/Documents/Databases/AILUN/"
historyFile <- "hpc111"
newPlatforms.list<-vector("list", 0)
for (i in 1:length(newPlatforms)){
   print(newPlatforms[i])
   conn <- gzfile(description = paste(newPlatforms[i], ".annot.gz", sep = ""), open = "rt")
   data <- read.delim(conn, stringsAsFactors = FALSE, header = FALSE)
   close(conn)
   annotation<-data.frame(EntrezID = data[,2], ID = data[,1], stringsAsFactors = FALSE)
   annotation$EntrezID<-entrezUpdate(annotation$EntrezID, historyFile = historyFile)
   # remove duplicated 
   annotation<-annotation[!duplicated(annotation$ID),]
   rownames(annotation)<-annotation$ID
   newPlatforms.list <- append(newPlatforms.list, list(list(ann = annotation)))
}
names(newPlatforms.list)<-newPlatforms

chipframe<-append(chipframe, newPlatforms.list)

save(chipframe, file = "/home/galtschu2/fingerprint/data/chipframe.RData")
   
# end of script