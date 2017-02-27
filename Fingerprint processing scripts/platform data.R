# Retrieve platform info for fingerprint chips
# Author: Gabriel Altschuler
# Timestamp: 20110502
# Updated: 20110810
# Status: Updated

# load chipframe data
# Use GEOmetaDB to retrieve the arrays

library(GEOmetadb)

# save GEOmetadb database locally
setwd("/Users/GabrielAltschuler/Documents/Databases/GEOmetadb")

# download and extract database
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

# load chipframe for fingerprint
#load("/Users/GabrielAltschuler/Dropbox/fingerprint/data/chipframe.RData")
library(pathprint.v0.3.beta4)
# Retrieve information from GEO
data(chipframe)
chipnames<-names(chipframe)
dbListFields(con, "gpl")

rs <- dbGetQuery(con, "select * from gpl")

chipnames.GEOdata<-rs[match(chipnames, rs$gpl),c("gpl", "manufacturer", "organism", "title")]
write.table(chipnames.GEOdata, file = "/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Manuscripts/Main publication/supplementary data/platforms.txt", sep = "\t", quote = FALSE, row.names = FALSE)
dbDisconnect(con)