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

fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/SCE/"
setwd(fingerpath)

# retrieve list of files
# 200 bytes used as cutoff for successful fingerprint calculation

files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
unprocessedFiles<-dir(fingerpath)[file.info(dir(fingerpath))$size < 200]
files<-files[grep("SCE_KEGG_Wiki_static_", files)]
unprocessedFiles<-unprocessedFiles[grep("SCE_KEGG_Wiki_static_", unprocessedFiles)]

length(files)
length(unprocessedFiles)

# 159343 successful
# 2903 unsucessful

# some unsuccessful files may be due to incomplete or interupted download from GEO

GEOpath <- "/home/galtschu2/Documents/Databases/GEOfiles/"
setwd(GEOpath)
GEOfiles<-dir(GEOpath)
GEOfiles.size<-file.info(dir(GEOpath))$size
head(sort(table(GEOfiles.size), decreasing = FALSE))

GEOfiles<-gsub(".soft", "", GEOfiles)
unprocessedFiles<-gsub("SCE_KEGG_Wiki_static_", "", unprocessedFiles)

names(GEOfiles.size)<-GEOfiles
head(table(GEOfiles.size[unprocessedFiles]),100)

# there are 328 GEO files with zero file size
# these files have now been removed from the GEO and fingerprint directory
# sourced GEOfilesdownload_server.R to re-aquire these GEO files and re-fingerprint

