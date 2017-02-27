# Retrieve GPL IDs and titles from the GEO database
# Author: Gabriel Altschuler
# Timestamp: 20110822
# Status: Complete

# Retrieve full list
library(GEOmetadb)
# set directory to hold database file
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")
# download database (only need to do this once)
getSQLiteFile()
# retrieve records
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
gpl.title<-dbGetQuery(con, "select gpl, title from gpl")
dbDisconnect(con)

# Cross reference these with platforms covered by pathprint

library(pathprint)
data(chipframe)
pathprint.gpl<-gpl.title[match(names(chipframe), gpl.title$gpl),]


###
# add names to chipframe - on local mac
library(pathprint.v0.3.beta4)
data(chipframe)

for (i in names(chipframe)){
	chipframe[[i]]$title <- gpl.title$title[match(i, gpl.title$gpl)]
		}

save(chipframe, file = "/Users/GabrielAltschuler/Dropbox/fingerprint/data/chipframe.RData")

# End