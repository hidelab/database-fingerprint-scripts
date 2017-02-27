# Addressing problem with the dimensions of fingerprint arrays
# Problem is that not all fingerprints are based on full expression arrays

# perhaps this data is available using GEOmetadb?

defineFile<-function(){
readline("enter filename: ")
}
dataPath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/"
dir(path = dataPath, pattern = "platform.frame.20")
print("select platform file")
platformfile<-defineFile()
load(paste(dataPath, platformfile, sep = ""))

# for pathprint v0.3 used platform.frame.2011-06-12.RData
# for pathprint v0.3beta3 used platform.frame.2011-06-22.RData


# download and extract metadata database
library(GEOmetadb)
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

# retrieve metadata
gsm.char = dbGetQuery(con, "select gsm, characteristics_ch1 from gsm")
gsm.title = dbGetQuery(con, "select gsm, title from gsm")
gsm.source = dbGetQuery(con, "select gsm, source_name_ch1 from gsm")
gsm.species = dbGetQuery(con, "select gsm, organism_ch1 from gsm")
gsm.nrow = dbGetQuery(con, "select gsm, data_row_count from gsm")
gsm.gse<-geoConvert(rownames(platform.frame), c("gse"))$gse
gsm.gpl<-geoConvert(rownames(platform.frame), c("gpl"))$gpl
gpl.title = dbGetQuery(con, "select gpl, title from gpl")
gpl.species<-dbGetQuery(con, "select gpl, organism from gpl")
dbDisconnect(con)

GEO.metadata.matrix<-as.data.frame(matrix(nrow = length(rownames(platform.frame)), ncol = 7))
colnames(GEO.metadata.matrix)<-c("GSM", "GSE", "GPL", "Species", "Title", "Source", "Characteristics")
GEO.metadata.matrix$GSM<-as.character(rownames(platform.frame))
GEO.metadata.matrix$GSE<-gsm.gse$to_acc[match(GEO.metadata.matrix$GSM, gsm.gse$from_acc)]
GEO.metadata.matrix$GPL<-gsm.gpl$to_acc[match(GEO.metadata.matrix$GSM, gsm.gpl$from_acc)]
GEO.metadata.matrix$Species<-gsm.species$organism_ch1[match(GEO.metadata.matrix$GSM, gsm.species$gsm)]
GEO.metadata.matrix$Title<-gsm.title$title[match(GEO.metadata.matrix$GSM, gsm.title$gsm)]
GEO.metadata.matrix$Source<-gsm.source$source_name_ch1[match(GEO.metadata.matrix$GSM, gsm.source$gsm)]
GEO.metadata.matrix$Characteristics<-gsm.char$characteristics_ch1[match(GEO.metadata.matrix$GSM, gsm.char$gsm)]
GEO.metadata.matrix$nrow<-gsm.nrow$data_row_count[match(GEO.metadata.matrix$GSM, gsm.nrow$gsm)]
GEO.metadata.matrix$platformSpecies<-gpl.species$organism[
  match(GEO.metadata.matrix$GPL, gpl.species$gpl)]


platforms<-levels(as.factor(GEO.metadata.matrix$GPL))
platforms.nrow<-lapply(platforms, function(x){table(GEO.metadata.matrix$nrow[GEO.metadata.matrix$GPL == x])})
names(platforms.nrow)<-platforms
maxProbes<-sapply(platforms.nrow, function(x){max(as.numeric(names(x)))})
# there are a large number of platforms that have less than the maximum
# is the maximum always the max
sapply(platforms.nrow, function(x){as.numeric(names(which.max(x)))}) == 
  sapply(platforms.nrow, function(x){max(as.numeric(names(x)))})

# apply a threshold defining that if an array ommits <10% of the probes it should not be fingerprinted 
GEO.metadata.matrix$probes<-NA
for (i in 1:length(platforms)){
  threshold<-0.9*maxProbes[i]
  GEO.metadata.matrix$probes[
    GEO.metadata.matrix$GPL == platforms[i]
    ] <- (GEO.metadata.matrix$nrow[
              GEO.metadata.matrix$GPL == platforms[i]
              ] >= threshold)
}
  
  
platforms.probes<-lapply(platforms, function(x){table(GEO.metadata.matrix$probes[GEO.metadata.matrix$GPL == x])})

# overwrite files
library(pathprint.v0.3.beta2)
data(GEO.fingerprint.matrix)
GEO.metadata.matrix<-GEO.metadata.matrix


GPL<-gpl.title[match(names(table(GEO.metadata.matrix$GPL)), gpl.title$gpl),]
rownames(GPL)<-GPL$gpl



# must remove Agilent arrays as these are two channel and cannot be analyzed in the same way
agilents<-GPL[grep("gilent", GPL$title),1]

# "GPL1708" "GPL4133" "GPL4134" "GPL6466" "GPL6480" "GPL6848" "GPL7202" "GPL887"  "GPL891"
# also remove B chips
bchips<-c("GPL340", "GPL97")

GEO.metadata.matrix$probes[GEO.metadata.matrix$GPL %in% agilents]<-FALSE
GEO.metadata.matrix$probes[GEO.metadata.matrix$GPL %in% bchips]<-FALSE
GEO.metadata.matrix.old<-GEO.metadata.matrix
GEO.metadata.matrix<-GEO.metadata.matrix[GEO.metadata.matrix$probes,]
# and remove chips where platformSpecies and species do not match
GEO.metadata.matrix<-GEO.metadata.matrix[
  GEO.metadata.matrix$Species == GEO.metadata.matrix$platformSpecies,]

GEO.metadata.matrix<-GEO.metadata.matrix[,-match("probes", colnames(GEO.metadata.matrix))]
rownames(GEO.metadata.matrix)<-1:nrow(GEO.metadata.matrix)
# now re-save
save(GEO.metadata.matrix, file = "/home/galtschu2/fingerprint/data/GEO.metadata.matrix.RData")
library(pathprint.v0.3.beta2)
data(GEO.fingerprint.matrix)
GEO.fingerprint.matrix<-GEO.fingerprint.matrix[,GEO.metadata.matrix$GSM]
save(GEO.fingerprint.matrix, file = "/home/galtschu2/fingerprint/data/GEO.fingerprint.matrix.RData")