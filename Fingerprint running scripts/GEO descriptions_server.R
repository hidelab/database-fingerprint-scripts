
GEOpath <- "/home/galtschu2/Documents/Databases/GEOfiles/"
GEOfiles <- gsub(".soft", "", dir(path = GEOpath))
	
# retrieve GEO file and strip out reuqired meta data
library(GEOquery)
description<-vector("list", 0)

for (i in 1:length(GEOfiles)){
	GSM<-GEOfiles[i]
	print("Downloading/loading GEO file")
	if (GSM %in% GEOfiles){
		try(geo <- getGEO(GSM, file = (paste(GEOpath, GSM, ".soft", sep = ""))))
		}
	else if (!(GSM %in% GEOfiles)){
		try(geo<-getGEO(GSM, destdir = GEOpath))
		}
	geoDescription<-list(Meta(geo)$description)
	names(geoDescription)<-GEOfiles[i]	
	description<-append(description, geoDescription)
	}

save(description, file = "/home/galtschu2/Documents/Databases/GEOdescription.RData")