# Fingerprinting worker script - running single chip enrichment option
# Delay allows time to monitor process starting up
print("preparing for script launch, 5s synch delay")
Sys.sleep(5)

# Read the genelist in
genelist <- sampleset[[7]];
total <- length(genelist)

# load contents of directory once
# this used to be within the loop
dir<-dir(path = path)

for (i in 1:length(genelist)){
	print(paste(i, "of", total, sep = " "))
  filename<-paste(header, genelist[i], sep = "_")
	if (filename %in% dir){
		print ("File already analyzed")
		  }
	else if (!(filename %in% dir(path = path))){		
		temp<-geo2Rank(
              GSM = genelist[i],
              GEOpath = GEOpath
              )
		temp1<-list("name" = temp)
		names(temp1)<-genelist[i]
		save(temp1, file = paste(path, filename, sep =""))

		  }
		}
