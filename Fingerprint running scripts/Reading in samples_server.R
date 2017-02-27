# Creating updated list of GEO references
#
# Author: Gabriel Altschuler
# Status: Updated
# Timestamp: 20111128
# Script to read list of GEO samples to update the fingerprint analysis with new GEO files
# The list of GEO files for each platform below is obtained from the GEO website
# GEOsamples is the final object that gets split up and passed to the fingerprinting scripts
# Currently split into 14 parts as the fingerprinting script is passed to 14 cores
# Each tackles a random subset of the arrays to prevent a bias for one particular platform
# This is required as some chips take longer to run than others

# load in chipfile if not already in global environment
# This contains information on the platforms to be included in the fingerprint
# Should already be in the global environment from preceeding script
if (!(exists("chipframe"))){
	running.dir <- dirname(parent.frame(2)$ofile)
	dataPath<- gsub("database-fingerprint-scripts/Fingerprint running scripts", "fingerprint/data", running.dir)
	if(running.dir == "."){
		dataPath <- "../../fingerprint/data"
		}
	load(paste(dataPath, "chipframe.RData", sep="/"), .GlobalEnv)
	}

# Use GEOmetaDB to retrieve the arrays
wd<-getwd()
library(GEOmetadb)

# save GEOmetadb database locally on hpc111
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")

# download and extract database
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

# Retrieve information from GEO
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

# reset working directory
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

