# Script to initiate an update of the pathway fingerprint
#
# Author: Gabriel Altschuler
# Status: Updated to allow source
# Timestamp: 20111202
# Location: hpc111
# Script to initiate an update of the pathway fingerprint 
# This script acts as a log of all of the fingerprint updates
# scroll to bottom for most recent and to add another


# BEFORE RUNNING THIS SCRIPT SEE IF CAN DELETE THE RAW FINGERPRINT FILES IN THE 
# RANKING DIRECTORY
# 2nd August 2013 - rank matrix - GPL570, GPL1261
setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
source("GEOsurvey_server_unionRank.R")


# 24th June 2013 - Union DPD genesets re-run using common genes only - GPL570, GPL1261
setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
source("GEOsurvey_server_unionDPD.R")

# 29th April 2013 - Union genesets re-run using common genes only - GPL570, GPL1261
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
#source("GEOsurvey_server_union.R")

# 8th April 2013 - Union genesets - GPL570, GPL1261
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
#source("GEOsurvey_server_union.R")

# 23rd December 2012 - Drug-Pathways-Diseases
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
#source("GEOsurvey_server_DPD.R")

# 10th April 2012 - running transcription factor modules
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint #running scripts/")
#source("GEOsurvey_server_TF.R")

# 2nd April 2012 - running mean rank version for Yered
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint #running scripts/")
#source("GEOsurvey_server_mean.R")

# 28th November - running squared rank version, adding Mouse and Human ST chips
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
#source("GEOsurvey_server_sq.R")

#
## 2nd August - running squared rank random genesets version
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
#source("GEOsurvey_server_sq_random.R")
#
## 15th July - running squared rank version
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
#source("GEOsurvey_server_sq.R")
#
#
## 11th June - running random genesets
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
#source("GEOsurvey_server_random.R")
#
## 8th June 2011 - adding model species and using updated genesets
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/")
#source("GEOsurvey_server_update.R")
#


##########
## running update script 4_10_10
#
## R server session
#
#setwd("/home/galtschu2/Documents/Projects/Fingerprinting/scripts")
#source("GEOfilesdownload_server.R")
#
## and in a new R server session
#
#setwd("/home/galtschu2/Documents/Projects/Fingerprinting/scripts")
#source("GEOsurvey_server.R")
#
## updating for the ABI chips
#
#setwd("/home/galtschu2/Documents/Projects/Fingerprinting/scripts")
#source("GEOsurvey_server_update.R")
#
## Running update for the experimental gene sets
#setwd("/home/galtschu2/Documents/Projects/Fingerprinting/scripts")
#source("GEOsurvey_server_experimental.R")
#
## this set of 10 pathways is running at a rate of 1 fingerprint every 2.2 seconds, will take approx 3 days to run out.
#
## Running ranking script - this runs at a rate of 1024 in 6 mins = (2.8 per second) therefore will take approx 12 hours to rank the full database
#setwd("/home/galtschu2/Documents/Projects/Ranking/scripts")
#source("GEOsurvey_server_ranking.R")
#
#
## implementing fingerprinting algorithm using single.chip.enrichment rather than single.chip.GSEA
## initially use same pathways and expression arrays to allow a direct comparison, therefore do not run the GEO download script
## need to edit source("GEOsurvey_server_update.R") to use new algorithm
## Don't forget to screen server terminal as will take a while to run!
#
## this was run on 23rd March 2010
#
#setwd("/home/galtschu2/Documents/Projects/Fingerprinting/scripts")
#source("GEOsurvey_server_update.R")
#
#
## Removed ~ 320 GEO files of zero length and re-run fingerprint
#setwd("/home/galtschu2/Documents/Projects/Fingerprinting/scripts")
#source("GEOsurvey_server_update.R")
#
## problem with one Agilent dataset, GSE12553
## the species data is entered into ch_2 rather that ch_1
## corrected script to error check for this
#
## now should be fine to run fingerprint again
## also updated to include additional HT Affy platforms
#
#setwd("/home/galtschu2/Documents/Projects/Fingerprinting/scripts")
#source("GEOsurvey_server_update.R")
#
## update 30th March 2010 has 19243 missing fingerprints - this seems quite a few!
## try to identify which node crashed
#for (i in 1:length(sampleset)){
#  missing <- sum(!(unlist(sampleset[[i]]) %in% dir))
#  print(paste("Node", i, "missing", missing, sep = " "))
#  if (missing > 0){
#    crashPoint <- head(sampleset[[i]][(!(unlist(sampleset[[i]]) %in% dir))],1)
#    print(paste("Node", i, "crash point", crashPoint, sep = " "))
#    }
#  }
#
#
## re-run to include the newly added Illumina and Agilent arrays 30th March 2010
## Before running need to update chipframe file to chipframe_update_29_3_10.RData
## updated worker scripts (SCE_parallel_a.R etc to only load directory once at the start
## this should speed this step up by about 1s per array
## Need to update
## 1)
## GSEA working script within "gabriel functions.R"
## sftp
## from /Users/GabrielAltschuler/Dropbox/fingerprint/scripts
## to /home/galtschu2/Documents/Databases
## 2)
## GEOsurvey_server_update.R
## Reading in samples_server_update.R
## sftp
## from /Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Scripts
## to /home/galtschu2/Documents/Projects/Fingerprinting/scripts 
## ssh to hpc111, screen terminal, start R and run
#
#setwd("/home/galtschu2/Documents/Projects/Fingerprinting/scripts")
#source("GEOsurvey_server_update.R")
#
#
## at the end of the run can assess which fingerprints were/were not processed
#path<-"data/Fingerprints/SCE/"
#dir<-dir(path)
#dir<-gsub("SCE_KEGG_Wiki_static_", "", dir)
#sum(!(unlist(sampleset) %in% dir))
#








