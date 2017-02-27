# Fingerprinting worker script - running single chip enrichment option
# Delay allows time to monitor process starting up
# Experimental genesets
print("preparing for script launch, 5s synch delay")
Sys.sleep(5)
# Read the genelist in
genelist <- sampleset[[10]];
total <- length(genelist)
# load contents of directory once
# this used to be within the loop
path<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/SCE/experimental/"
dir<-dir(path = path)
# removed line to record vector as a waste of memory and not used
# vector<-vector("list", 0)

for (i in 1:length(genelist)){
	print(paste(i, "of", total, sep = " "))
  filename<-paste("SCE_Experimental", genelist[i], sep = "_")
	if (filename %in% dir){
		print ("File already analyzed")
		  }
	else if (!(filename %in% dir(path = path))){		
		temp<-geo2fingerprint(
              GSM = genelist[i], 
              GEOthreshold = FALSE,
              GEOpath = "/home/galtschu2/Documents/Databases/GEOfiles/",
              geneset = "Experimental",
              enrichmentMethod = "SCE",
              transformation = "rank",
              statistic = "mean",
              normalizedScore = FALSE,
              progressBar = FALSE
              )
		temp1<-list("name" = temp)
		names(temp1)<-genelist[i]
		save(temp1, file = paste(path, filename, sep =""))
#   vector<-append(vector, temp1)
		  }
		}
# save(vector, file="data/parallel_vector1.R")
