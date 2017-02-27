# downloading GEOfiles for fingerprinting
# at the moment this runs as a stand-alone script
# should be integrated into the fingerprinting process
# The main fingerprint script geo2fingerprint has been update
# The directory of GEO files is now only searched if the directory list is now in the global environment
# If this script is to be run concurrently then need to insert the directory list into the same global environment

source("/home/galtschu2/Documents/Projects/Fingerprinting/scripts/Reading in samples_server_update.R")
library(GEOquery)
GEOpath <- "/home/galtschu2/Documents/Databases/GEOfiles/"
GEOfiles <- gsub(".soft", "", dir(path = GEOpath))

genelist <- unlist(sampleset);
	for (i in 1:length(genelist)){
		GSM<-genelist[i]
		if (GSM %in% GEOfiles){
			print("file already downloaded")
			}
		else if (!(GSM %in% GEOfiles)){
			try(temp<-getGEO(GSM, destdir = GEOpath))
			}
		}
