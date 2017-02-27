# Script to read list of GEO samples to update the fingerprint analysis with new GEO files
# The list of GEO files for each platform below is obtained from the GEO website
# GEOsamples.new is the final object that gets passed to the fingerprinting scripts
# This script contains sequential updates, added in series, most recent from March 29th 2010

print("Compiling list of GEO files to fingerprint")

setwd("/home/galtschu2/Documents/Databases/PlatformFiles/10_7_10/")

GPL570.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL570.txt")[grep("!Platform_sample_id", readLines("GPL570.txt"))])
GPL1261.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL1261.txt")[grep("!Platform_sample_id", readLines("GPL1261.txt"))])
GPL339.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL339.txt")[grep("!Platform_sample_id", readLines("GPL339.txt"))])
GPL96.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL96.txt")[grep("!Platform_sample_id", readLines("GPL96.txt"))])
GPL97.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL97.txt")[grep("!Platform_sample_id", readLines("GPL97.txt"))])
GPL81.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL81.txt")[grep("!Platform_sample_id", readLines("GPL81.txt"))])
GPL8321.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL8321.txt")[grep("!Platform_sample_id", readLines("GPL8321.txt"))])
GPL8300.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL8300.txt")[grep("!Platform_sample_id", readLines("GPL8300.txt"))])
GPL340.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL340.txt")[grep("!Platform_sample_id", readLines("GPL340.txt"))])
GPL571.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL571.txt")[grep("!Platform_sample_id", readLines("GPL571.txt"))])

GEOsamples.10_7_10<-c(GPL570.samples, GPL1261.samples, GPL339.samples, GPL96.samples, GPL97.samples, GPL81.samples, GPL8321.samples, GPL8300.samples, GPL340.samples, GPL571.samples)

setwd("/home/galtschu2/Documents/Databases/PlatformFiles/4_10_10/")

GPL570.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL570.txt")[grep("!Platform_sample_id", readLines("GPL570.txt"))])
GPL1261.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL1261.txt")[grep("!Platform_sample_id", readLines("GPL1261.txt"))])
GPL339.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL339.txt")[grep("!Platform_sample_id", readLines("GPL339.txt"))])
GPL96.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL96.txt")[grep("!Platform_sample_id", readLines("GPL96.txt"))])
GPL97.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL97.txt")[grep("!Platform_sample_id", readLines("GPL97.txt"))])
GPL81.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL81.txt")[grep("!Platform_sample_id", readLines("GPL81.txt"))])
GPL8321.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL8321.txt")[grep("!Platform_sample_id", readLines("GPL8321.txt"))])
GPL8300.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL8300.txt")[grep("!Platform_sample_id", readLines("GPL8300.txt"))])
GPL340.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL340.txt")[grep("!Platform_sample_id", readLines("GPL340.txt"))])
GPL571.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL571.txt")[grep("!Platform_sample_id", readLines("GPL571.txt"))])
GPL2986.samples<-gsub("!Platform_sample_id = ", "", readLines("GPL2986.txt")[grep("!Platform_sample_id", readLines("GPL2986.txt"))])

GEOsamples.4_10_10<-c(GPL570.samples, GPL1261.samples, GPL339.samples, GPL96.samples, GPL97.samples, GPL81.samples, GPL8321.samples, GPL8300.samples, GPL340.samples, GPL571.samples)

# In the update 10_7_10 only wanted to add the new fingerprints
# 4489 new samples
GEOsamples.new<-setdiff(GEOsamples.4_10_10, GEOsamples.10_7_10)

# In the ABI array update only wanted to fingerprint the ABI arrays
# fingerprinting ABI array
GEOsamples.new<-GPL2986.samples

# In the gene.set.enrichment update (23rd March 2011) want to fingerprint all arrays
GEOsamples.new<-union(GPL2986.samples, GEOsamples.10_7_10)

# further update would be to use GEOmetaDB to retrieve the arrays
wd<-getwd()
library(GEOmetadb)
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")
# download and extract database
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

# test<-dbGetQuery(con, "SELECT gsm, gpl FROM gsm WHERE gpl = 'GPL2986'")
# test<-dbGetQuery(con, "SELECT gsm, gpl FROM gsm WHERE gpl = 'GPL2986'")

# load in chipfile
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe.RData")
chipnames<-names(chipframe)

chipvector<-vector("character", 0)
print("Querying GEOmetaDB")
for (i in 1:length(chipnames)){
  query<-paste("SELECT gsm FROM gsm WHERE gpl = '", chipnames[i], "'", sep = "")
  print(query)
  gsmList<-dbGetQuery(con, query)
  chipvector<-c(chipvector, unlist(gsmList, use.names = FALSE))
  }
  
GEOsamples.new<-chipvector

dbDisconnect(con)
setwd(wd)

# scramble sample names to avoid bias in processing time asscociated with the larger chips
# apportion into 12 or 14
print(paste(length(GEOsamples.new), "to be fingerprinted", sep = " ")) 
GEOsamples.rand <- sample(GEOsamples.new)
len <- length(GEOsamples.rand) %/% 14

sampleset<-vector("list", 14)
for (i in 1:length(sampleset)){
	sampleset[[i]] <- GEOsamples.rand[(1+((i-1)*len)):(i*len)]
	}
sampleset[[14]] <- GEOsamples.rand[(1+(13*len)):length(GEOsamples.rand)]

