# Post processing fingerprint_step 1
# Author: Gabriel Altschuler
# Timestamp: 20130610
# Status: Re-compiled to be sourced
# location: hpc111
#
# Processing the fingerprinted data
# Probability of expression method
# following Parmigiani et al 2002
# based on "analysis of SCG histograms_server.R"

# SCRIPT IS IN 3 DISTINCT PARTS, THIS IS 1 OF 3
# quit and restart R between each to conserve memory
# It is not properly returned to system otherwise (even with garbage collection)

########
# PART 1
########

# First require matrix of raw single chip enrichments
# Using SCE processed data, mean rank method, DPD genesets

# set directory holding fingerprint collection
fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/UnionDPD/"
# define file identifier
header<-"UnionDPD_"
setwd(fingerpath)

# retrieve list of files, some will not have been properly processed
# use 200 bytes used as cutoff for successful fingerprint calculation

print("Loading files")

files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
unprocessedFiles<-dir(fingerpath)[file.info(dir(fingerpath))$size < 200]
files<-files[grep(header, files)]
unprocessedFiles<-unprocessedFiles[grep(header, unprocessedFiles)]
print(paste(length(files), "successfully fingerprinted GEO arrays"))
print(paste(length(unprocessedFiles), "GEO arrays could not be fingerprinted"))


############
# Remove arrays with <90% probe representation
print("Removing arrays with <90% probe representation")
library(GEOmetadb)
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")
#getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")

# retrieve metadata, including number of rows in the GEO data table
gsm.nrow = dbGetQuery(con, "select gsm, data_row_count from gsm")
gsm.gpl<-geoConvert(gsub(header, "", files), c("gpl"))$gpl
dbDisconnect(con)

probe.data<-as.data.frame(matrix(nrow = length(files), ncol = 3))
colnames(probe.data)<-c("GSM","GPL","nrow")
probe.data$GSM<-gsub(header, "", files)
probe.data$GPL<-gsm.gpl$to_acc[match(probe.data$GSM, gsm.gpl$from_acc)]
probe.data$nrow<-gsm.nrow$data_row_count[match(probe.data$GSM, gsm.nrow$gsm)]

# for each platform build a table of the data table nrows and find the max (i.e. complete array)
platforms<-levels(as.factor(probe.data$GPL))
platforms.nrow<-lapply(platforms, function(x){table(probe.data$nrow[probe.data$GPL == x])})
names(platforms.nrow)<-platforms
maxProbes<-sapply(platforms.nrow, function(x){max(as.numeric(names(x)))})

# define a TRUE/FALSE column for whether an array has at least 90% of the probes accounted for
probe.data$probes<-FALSE
probe.data$GPL[is.na(probe.data$GPL)]<- "Missing"
for (i in 1:length(platforms)){
  threshold<-0.9*maxProbes[i]
  probe.data$probes[
    probe.data$GPL == platforms[i]
    ] <- (probe.data$nrow[
              probe.data$GPL == platforms[i]
              ] >= threshold)
}

# in addition, remove the Agilent arrays
# must remove Agilent arrays as these are two channel and cannot be analyzed in the same way
agilents<-c("GPL1708", "GPL4133", "GPL4134", "GPL6466", "GPL6480", "GPL6848", "GPL7202", "GPL887",  "GPL891")
probe.data$probes[probe.data$GPL %in% agilents]<-FALSE

# now re-define the file list and continue
files<-paste(header, probe.data$GSM[probe.data$probes], sep = "")
setwd(fingerpath)
length(files)
#######################
# HUMAN CORRELATION MATRIX
humanFiles <- paste(header, probe.data$GSM[
	probe.data$probes & probe.data$GPL == "GPL570"
											], sep = "")
load(humanFiles[1])
pathways <- rownames(temp1[[1]]$SCG)
# NEED TO SPLIT INTO HUMAN AND MOUSE AND PROCESS SEPARATELY
#load("validPathwayIndex.R")
print("Creating human expression matrix")
# pre-allocate memory for matrix
SCE.matrix <- matrix(nrow = length(humanFiles), ncol = length(pathways))
#SCE.matrix <- matrix(nrow = length(humanFiles), ncol = sum(validPathwayIndex))
colnames(SCE.matrix)<-pathways
rownames(SCE.matrix)<-gsub(header, "", humanFiles)
print(object.size(SCE.matrix), units = "auto")
SCE.matrix[,] <- as.integer(1)
pb <- txtProgressBar(min = 0, max = length(humanFiles), style = 3)
for (i in 1:length(humanFiles)){
  load(humanFiles[i])
  print(i)
  # setTxtProgressBar(pb, i)
  tempData <- as.integer(temp1[[1]]$SCG)
  #tempData <- tempData[validPathwayIndex]
  tempData[is.na(tempData)] <- 0
  SCE.matrix[i,]<-as.integer(tempData)
  if (i %% 1000 == 0){
  	gc()
  	}
  }

# remove cols that are all NAs
# straight forward method does not work
invalidCountCol <- colSums(SCE.matrix == 0)
SCE.matrix <- SCE.matrix[,invalidCountCol < nrow(SCE.matrix)]
invalidCountRow <- rowSums(SCE.matrix == 0)
# remove all rows that are not complete - probably something wrong with the way that data that has been uploaded
SCE.matrix <- SCE.matrix[invalidCountRow == 0,]

# check that all ok
sum(SCE.matrix == 0)

# save matrix
print("Saving matrix. The file size is")
print(object.size(SCE.matrix), units = "auto")
save(SCE.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, "GPL570.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

######################
# calculate human correlation matrix
# create line by line
# multiply by 10000000 to maintain integer values

cor.matrix <- matrix(nrow = ncol(SCE.matrix),
					 ncol = ncol(SCE.matrix))
colnames(cor.matrix)<-colnames(SCE.matrix)
rownames(cor.matrix)<-colnames(SCE.matrix)
print(object.size(cor.matrix), units = "auto")


imax = floor(ncol(SCE.matrix)/20) + 1		 
for (i in 1:imax){
  print(i)
  istart = 1 + 20 * (i-1)
  iend = istart + 19
  if (iend > ncol(SCE.matrix)){
  	iend = ncol(SCE.matrix)
  	}
  # setTxtProgressBar(pb, i)
  cor.matrix[,istart:iend]<-as.integer(1000000*cor(
  						SCE.matrix,
  						SCE.matrix[,istart:iend]),
  						use = "all.obs")
  						
  	gc()
  }		


save(cor.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, ".human.cor.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

############
#######################
# MOUSE CORRELATION MATRIX
mouseFiles <- paste(header, probe.data$GSM[
	probe.data$probes & probe.data$GPL == "GPL1261"
											], sep = "")
load(mouseFiles[1])
pathways <- rownames(temp1[[1]]$SCG)
print("Creating mouse expression matrix")
# pre-allocate memory for matrix
SCE.matrix <- matrix(nrow = length(mouseFiles), ncol = length(pathways))
colnames(SCE.matrix)<-pathways
rownames(SCE.matrix)<-gsub(header, "", mouseFiles)
print(object.size(SCE.matrix), units = "auto")
SCE.matrix[,] <- as.integer(1)
pb <- txtProgressBar(min = 0, max = length(mouseFiles), style = 3)
for (i in 1:length(mouseFiles)){
  load(mouseFiles[i])
  print(i)
  # setTxtProgressBar(pb, i)
  tempData <- as.integer(temp1[[1]]$SCG)
  #tempData <- tempData[validPathwayIndex]
  tempData[is.na(tempData)] <- 0
  SCE.matrix[i,]<-as.integer(tempData)
  if (i %% 100 == 0){
  	gc()
  	}
  }

# remove cols that are all NAs
# straight forward method does not work
invalidCountCol <- colSums(SCE.matrix == 0)
SCE.matrix <- SCE.matrix[,invalidCountCol < nrow(SCE.matrix)]
invalidCountRow <- rowSums(SCE.matrix == 0)
# remove all rows that are not complete - probably something wrong with the way that data that has been uploaded
SCE.matrix <- SCE.matrix[invalidCountRow == 0,]

# check that all ok
sum(SCE.matrix == 0)

# save matrix
print("Saving matrix. The file size is")
print(object.size(SCE.matrix), units = "auto")
save(SCE.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, "GPL1261.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )



######################
# calculate mouse correlation matrix
# create line by line
# multiply by 10000000 to maintain integer values

cor.matrix <- matrix(nrow = ncol(SCE.matrix),
					 ncol = ncol(SCE.matrix))
colnames(cor.matrix)<-colnames(SCE.matrix)
rownames(cor.matrix)<-colnames(SCE.matrix)
print(object.size(cor.matrix), units = "auto")


imax = floor(ncol(SCE.matrix)/20) + 1		 
for (i in 1:imax){
  print(i)
  istart = 1 + 20 * (i-1)
  iend = istart + 19
  if (iend > ncol(SCE.matrix)){
  	iend = ncol(SCE.matrix)
  	}
  # setTxtProgressBar(pb, i)
  cor.matrix[,istart:iend]<-as.integer(1000000*cor(
  						SCE.matrix,
  						SCE.matrix[,istart:iend]),
  						use = "all.obs")
  						
  	gc()
  }		


save(cor.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, ".mouse.cor.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

# END OF SCRIPT
################


# test two human matrices
#GPL570a.cor.matrix <- cor.matrix
#GPL570b.cor.matrix <- cor.matrix
#GPL570.cor.matrix <- cor.matrix
#GPL1261.cor.matrix <- cor.matrix

#SCE.matrix.GPL570 <- SCE.matrix[
#	rownames(SCE.matrix) %in% probe.data$GSM[probe.data$GPL == "GPL570"]
#			,]


#SCE.matrix.GPL1261 <- SCE.matrix[
#	rownames(SCE.matrix) %in% probe.data$GSM[probe.data$GPL == "GPL1261"]
#			,]


#rm(cor.matrix)

imax = floor(ncol(SCE.matrix)/20) + 1		 
for (i in 1:imax){
  print(i)
  istart = 1 + 20 * (i-1)
  iend = istart + 19
  if (iend > ncol(SCE.matrix)){
  	iend = ncol(SCE.matrix)
  	}
  # setTxtProgressBar(pb, i)
  GPL570.cor.matrix[,istart:iend]<-as.integer(1000000*cor(
  						SCE.matrix.GPL570,
  						SCE.matrix.GPL570[,istart:iend]),
  						use = "all.obs")
  						
  GPL1261.cor.matrix[,istart:iend]<-as.integer(1000000*cor(
  						SCE.matrix.GPL1261,
  						SCE.matrix.GPL1261[,istart:iend]),
  						use = "all.obs")
  # split GPL matrix into two
  GPL570a.cor.matrix[,istart:iend]<-as.integer(1000000*cor(
  						SCE.matrix.GPL570[1:floor(nrow(SCE.matrix.GPL570)/2),],
  						SCE.matrix.GPL570[1:floor(nrow(SCE.matrix.GPL570)/2),istart:iend]),
  						use = "all.obs")
  GPL570b.cor.matrix[,istart:iend]<-as.integer(1000000*cor(
  						SCE.matrix.GPL570[(1+floor(nrow(SCE.matrix.GPL570)/2)):nrow(SCE.matrix.GPL570),],
  						SCE.matrix.GPL570[(1+floor(nrow(SCE.matrix.GPL570)/2)):nrow(SCE.matrix.GPL570),istart:iend]),
  						use = "all.obs")						

  	gc()
  }					 

save(cor.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, ".cor.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

save(GPL570.cor.matrix,
	 file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/",
	 			  header, ".GPL570.cor.matrix.",
     			  Sys.Date(), ".RData",
     			  sep = ""
     )
   )

save(GPL1261.cor.matrix,
	 file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/",
	 			  header, ".GPL1261.cor.matrix.",
     			  Sys.Date(), ".RData",
     			  sep = ""
     )
   )


# produce matrix of correlation differences

corr.diff <- (GPL570.cor.matrix - GPL1261.cor.matrix)/1000000

corr.diff.abs <- (GPL570.cor.matrix - GPL1261.cor.matrix)/abs(GPL570.cor.matrix)

which(abs(corr.diff[!grepl("Netpath", colnames(corr.diff)),
				!grepl("Netpath", colnames(corr.diff))
				]) > 0.5 &
				abs(corr.diff.abs[!grepl("Netpath", colnames(corr.diff)),
				!grepl("Netpath", colnames(corr.diff))
				]) > 2, arr.ind = TRUE)

# diseases are [181:227]
signif <- which(abs(corr.diff[!grepl("Netpath", colnames(corr.diff)),
				!grepl("Netpath", colnames(corr.diff))
				]) > 0.5 &
				abs(corr.diff.abs[!grepl("Netpath", colnames(corr.diff)),
				!grepl("Netpath", colnames(corr.diff))
				]) > 0.5, arr.ind = TRUE)				

signif <- which(abs(corr.diff) > 0.5 &
				abs(corr.diff.abs) > 0.5, arr.ind = TRUE)				
signifDisease <- signif[signif[,1] %in% 181:227,]
signifDisease <- cbind(signifDisease, colnames(corr.diff)[signifDisease[,2]])		

signifDisease <- cbind(signifDisease,
	apply(signifDisease, 1, function(x){GPL570.cor.matrix[as.numeric(x[1]), as.numeric(x[2])]/1000000}
	))
signifDisease <- cbind(signifDisease,
	apply(signifDisease, 1, function(x){GPL1261.cor.matrix[as.numeric(x[1]), as.numeric(x[2])]/1000000}
	))
colnames(signifDisease)[3:5] <- c("Pathway2", "HumanCor", "MouseCor")
	
	
write.table(signifDisease, file = "Results/signifDisease.txt", row.names = T, col.names = T, sep ="\t", quote = F)	
		
corr.diff.test <- (GPL570a.cor.matrix - GPL570b.cor.matrix)/1000000
p1 <- hist(corr.diff, 100)
p2 <- hist(corr.diff.test, 100)
plot( p1, col=rgb(0,0,1,1/4), xlim=c(-1.2,1.2))  # first histogram
plot( p2, col=rgb(1,0,0,1/4), xlim=c(-1.2,1.2), add=T)  # second
par(mfcol = c(2,1))
geo.pluripotentDistance.hist<-hist(geo.pluripotentDistance[,"distance"],
	nclass = 50, xlim = c(0,1), main = "Distance from pluripotent consensus")
par(mar = c(7, 4, 4, 2))
hist(geo.pluripotentDistance[pluripotents.frame$GSM, "distance"],
	breaks = geo.pluripotentDistance.hist$breaks, xlim = c(0,1), 
	main = "", xlab = "above: all GEO, below: curated pluripotent samples")

hist(corr.diff.test,100)