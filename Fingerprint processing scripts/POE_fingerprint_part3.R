# Normalizing fingerprint using POE_part3
#
# Author: Gabriel Altschuler
# Timestamp: 20110405
#
# to be run on hpc111
#
# Processing the fingerprinted data
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# This script had to be run in 3 parts to restore memory properly to the system

# each pathawy has been separately normalized for each frame
# now need to compile all data together into one large dataframe

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/")

# load list of platforms
load(grep("platform", dir(), value = TRUE))
platforms<-as.character(unique(platform.frame$Platform))

# load list of pathways
load("/home/galtschu2/fingerprint/data/kegg_wiki_TR_static.Hs.gs")
geneset<-kegg_wiki_TR_static.Hs.gs

# setup POE.matrix
POE.matrix<-matrix(nrow = length(geneset), ncol = length(platform.frame$Platform))
rownames(POE.matrix)<-names(geneset)
colnames(POE.matrix)<-gsub("SCE_", "", platform.frame$GEO)

pb <- txtProgressBar(min = 0, max = length(geneset), style = 3)
for (j in 1:length(platforms)){
  array<-platforms[j]
  print(paste("Inserting POE values for ", array, sep = ""))
  for (i in 1:length(geneset)){
    load(paste("POE/pathway.SCG", i, array, "RData", sep = "."))
    POE.matrix[i,gsub("SCE_", "", names(pathway.POE))]<-pathway.POE
    setTxtProgressBar(pb, i)
    }
  }

save(POE.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/POE.matrix",
                              Sys.Date(), "RData",
                              sep = "."
                              )
    )

    
# save matrix with threshold applied at 0.5
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")
threshold<-seq(0.1,0.9,0.1)
#for (i in 1:length(threshold)){
for (i in 5){
  high<-threshold[i]
	low<-(-threshold[i])
	POE.matrix.threshold<-((POE.matrix)>high)-((POE.matrix)<low)
	POE.matrix.threshold[is.na(POE.matrix.threshold)]<-0
	assign(paste("POE.matrix.threshold.", threshold[i], sep = ""), POE.matrix.threshold)
	save(	list = paste("POE.matrix.threshold.", threshold[i], sep = ""),
			file = paste("POE.matrix.threshold.", threshold[i], ".RData", sep = "")
			)
	}

####################################

# load dataframe and threshold


    
  
  



