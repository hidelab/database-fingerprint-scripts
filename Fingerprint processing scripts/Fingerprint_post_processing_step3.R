# Post processing fingerprint_step 3
# Author: Gabriel Altschuler
# Timestamp: 20111130
# Status: Re-compiled to be sourced
# location: hpc111
#
# Processing the fingerprinted data
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# SCRIPT IS IN 3 DISTINCT PARTS, THIS IS 3 OF 3
# quit and restart R between each to conserve memory
# It is not properly returned to system otherwise (even with garbage collection)
########
# PART 3
########


# define paths and headers
header<-"sq_"
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")

Pathwaypath <- "/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/"
POEpath <- "/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/POE/"

# Define the location of the pathway frame
# If this script is run on same day as part 1  the line below loads the correct file
# If day is different, find the right files to load and adjust accordingly
# load(paste("platform.frame", Sys.Date(), "RData", sep = "."))
dir()[grep("platform", dir())]
filename <- readline(prompt = "Enter platform filename ") 
print("Loading platform data")
load(filename)

##########
# Convert platform names
platform.frame$Platform<-gsub("GPL96", "HGU133A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL571", "HGU133A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL3921", "HGU133A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL4685", "HGU133A", platform.frame$Platform)

platform.frame$Platform<-gsub("GPL4133", "G4112", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL6480", "G4112", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL1708", "G4112", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL6848", "G4112", platform.frame$Platform)

platform.frame$Platform<-gsub("GPL891", "G4121A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL6466", "G4121A", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL4134", "G4122F", platform.frame$Platform)
platform.frame$Platform<-gsub("GPL7202", "G4122F", platform.frame$Platform)

table(platform.frame$Platform)
platforms<-as.character(unique(platform.frame$Platform))
##########

# load list of pathways
library(pathprint)
pathwayNames<-names(get(genesets[[1]]))

# setup POE.matrix
POE.matrix<-matrix(nrow = length(pathwayNames), ncol = length(platform.frame$Platform))
rownames(POE.matrix)<-pathwayNames
colnames(POE.matrix)<-platform.frame$GEO

pb <- txtProgressBar(min = 0, max = length(pathwayNames), style = 3)
for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("Inserting POE values for ", array, sep = ""))
  for (i in 1:length(pathwayNames)){
    load(paste(POEpath, header, "/pathway.SCE.", i, ".", array, ".RData", sep = ""))
    POE.matrix[i,names(pathway.POE)]<-pathway.POE
    setTxtProgressBar(pb, i)
    }
  }

# remove arrays that are NA for the first pathway - this is usually due to a mis-match between the platform species and the annotated species
POE.matrix<-POE.matrix[,!(is.na(POE.matrix[1,]))]

save(POE.matrix, file = paste(header, ".POE.matrix.", Sys.Date(), ".RData", sep = ""))


POE.matrix.0.001<-(POE.matrix > 0.001) - (POE.matrix < -0.001)
# how many NAs? - These are probably due to the species or chips where a pathway is not represented properly
(sum(is.na(POE.matrix.0.001)))/(dim(POE.matrix.0.001)[1]*dim(POE.matrix.0.001)[2])
# how many zeros?
(sum((POE.matrix.0.001 == 0), na.rm = TRUE))/(dim(POE.matrix.0.001)[1]*dim(POE.matrix.0.001)[2])

save(POE.matrix.0.001, file = paste("/data/shared/Fingerprint/", header, ".POE.matrix.0.001.",
                              Sys.Date(), ".RData",
                              sep = ""
                              )
    )

#################
