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
fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/UnionRank/"
# define file identifier
header<-"UnionRank_"
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
getSQLiteFile()
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
humanFiles <- paste(header, probe.data$GSM[
	probe.data$probes & probe.data$GPL == "GPL570"
											], sep = "")
load(humanFiles[1])
pathways <- rownames(temp1[[1]]$SCG)
print("Creating human gene expression matrix")
# pre-allocate memory for matrix
# NEED TO MAKE SURE HERE THAT THE GENE NAMES ARE THE SAME FOR HUMAN AND MOUSE
# different gene names for human and mouse but this is ok as will be processed separately
# therefore need to re-order
# ALSO NEED TO GO BACK AND DELETE ALL OF THE FILES AFTERWARDS
SCE.matrix <- matrix(nrow = length(humanFiles), ncol = length(pathways))
#SCE.matrix <- matrix(nrow = length(humanFiles), ncol = sum(validPathwayIndex))
colnames(SCE.matrix)<-pathways
rownames(SCE.matrix)<-gsub(header, "", humanFiles)
SCE.matrix[,] <- as.integer(1)
print(object.size(SCE.matrix), units = "auto")
pb <- txtProgressBar(min = 0, max = length(humanFiles), style = 3)
for (i in 1:length(humanFiles)){
  load(humanFiles[i])
  # setTxtProgressBar(pb, i)
 # tempData <- rep(0,ncol(SCE.matrix))
 # if (all.equal(colnames(SCE.matrix), rownames(temp1[[1]]$SCG))){
 	#tempVector <- as.vector(temp1[[1]]$SCG)
  	tempData <- as.integer(temp1[[1]]$SCG)
  #tempData <- tempData[validPathwayIndex]
 	tempData[is.na(tempData)] <- 0
 	#tempData <- tempData
  #	}
  try(SCE.matrix[i,]<-as.integer(tempData))
  if (i %% 1000 == 0){
  	print(i)
  	gc()
  	}
  }
 
 invalidRow <- rowSums(SCE.matrix == 1) == ncol(SCE.matrix)
 
 SCE.matrix <- SCE.matrix[!invalidRow,]
 
# CHECK THAT INTEGER MATRIX BEFORE SAVING
# ALSO CHECK FILE SIZE
# IS THERE ANY WAY OF MAKING THIS FILE SMALLER?
# PERHAPS AS PRINCIPAL COMPONENTS BUT I DON'T THINK THIS WOULD PROVIDE THE NECESSARY ACCURACY
# COULD CHECK A FEW GENESETS HERE TO SEE IF GET THE SAME VALUES
# UP TO HERE
# NEXT STEP IS TO GET THIS INTO A PACKAGE AND TRY IT WITH SHINY
# 400 Mb?
# HOW MUCH RAM IS 
print("Saving matrix. The file size is")
print(object.size(SCE.matrix), units = "auto")
save(SCE.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, "GPL570.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

## UP TO HERE
# NEED TO
# 1) Calculate an example geneset for an example GSM and see if we get the same answer
# 2) Check that matrix is integer, is there a way to reduce the size further?
# perhaps need to calculate the data in parts?
######
# TESTING
single.chip.enrichment.preRanked <- function (exprs, geneset)
# single chip enrichment based on median rank
# pre-ranked data 
{
    Ns <- ncol(exprs)
    gene.names <- rownames(exprs)
    geneset.names <- names(geneset)
    score.matrix <- matrix(0, nrow = length(geneset), ncol = Ns)
    for (i in 1:length(geneset)) {
        overlap <- intersect(geneset[[i]], gene.names)
        if (length(overlap) == 0) {
            score.matrix[i, ] <- NA
        }
        else {
              score.matrix[i, ] <- apply(exprs, 2, function(x) {
                median(x[overlap])
                })
    		}
    	}
    colnames(score.matrix) <- colnames(exprs)
    rownames(score.matrix) <- geneset.names
    return(score.matrix)
}

setwd("~/database-fingerprint-scripts/unionGenesets")
# load("DPDUnion.genesets.RData")
# lapply(genesets, load, envir = .GlobalEnv)
# load chipframe with alternative GPL570 and GPL1261 annotations using only the commom genes
load("~/database-fingerprint-scripts/unionGenesets/chipframeCommon_GPL570_GPL1261.RData")
load("~/database-fingerprint-scripts/unionGenesets/DPD.Union.Hs.gs")
# calculate 1st geneset in first 5 arrays
testEnrichment <- single.chip.enrichment.preRanked(t(SCE.matrix[1:5,]),DPD.Union.Hs.gs)
# compare to values calcualted directly
testEnrichmentDirect <- matrix(nrow = length(DPD.Union.Hs.gs), ncol = 5)
rownames(testEnrichmentDirect) <- names(DPD.Union.Hs.gs)
colnames(testEnrichmentDirect) <- rownames(SCE.matrix)[1:5]
for (i in 1:5){
	tempfile = paste("~/Documents/Projects/Fingerprinting/data/Fingerprints/UnionDPD/UnionDPD_",
					 rownames(SCE.matrix)[i],
					 sep = "")
	load(tempfile)
	testEnrichmentDirect[,i]<-temp1[[1]]$SCG
	}

head(testEnrichment)	
head(testEnrichmentDirect)
cor(testEnrichment, testEnrichmentDirect,use = "complete.obs")
# looks like this has worked!! :)
# While 'here' calculate Jiantao's signature


#######
# Mouse signatures
#######################
mouseFiles <- paste(header, probe.data$GSM[
	probe.data$probes & probe.data$GPL == "GPL1261"
											], sep = "")
load(mouseFiles[1])
pathways <- rownames(temp1[[1]]$SCG)
print("Creating mouse gene expression matrix")
# pre-allocate memory for matrix
# NEED TO MAKE SURE HERE THAT THE GENE NAMES ARE THE SAME FOR HUMAN AND MOUSE
# different gene names for human and mouse but this is ok as will be processed separately
# therefore need to re-order
# ALSO NEED TO GO BACK AND DELETE ALL OF THE FILES AFTERWARDS
SCE.matrix <- matrix(nrow = length(mouseFiles), ncol = length(pathways))
colnames(SCE.matrix)<-pathways
rownames(SCE.matrix)<-gsub(header, "", mouseFiles)
SCE.matrix[,] <- as.integer(1)
print(object.size(SCE.matrix), units = "auto")
pb <- txtProgressBar(min = 0, max = length(mouseFiles), style = 3)
for (i in 1:length(mouseFiles)){
  load(mouseFiles[i])
  	tempData <- as.integer(temp1[[1]]$SCG)
 	tempData[is.na(tempData)] <- 0
  try(SCE.matrix[i,]<-as.integer(tempData))
  if (i %% 1000 == 0){
  	print(i)
  	gc()
  	}
  }
 
invalidRow <- rowSums(SCE.matrix == 1) == ncol(SCE.matrix)
 
SCE.matrix <- SCE.matrix[!invalidRow,]
 
# CHECK THAT INTEGER MATRIX BEFORE SAVING
# ALSO CHECK FILE SIZE
# IS THERE ANY WAY OF MAKING THIS FILE SMALLER?
# PERHAPS AS PRINCIPAL COMPONENTS BUT I DON'T THINK THIS WOULD PROVIDE THE NECESSARY ACCURACY
# COULD CHECK A FEW GENESETS HERE TO SEE IF GET THE SAME VALUES
# UP TO HERE
# NEXT STEP IS TO GET THIS INTO A PACKAGE AND TRY IT WITH SHINY
# 400 Mb?
# HOW MUCH RAM IS 
print("Saving matrix. The file size is")
print(object.size(SCE.matrix), units = "auto")
save(SCE.matrix, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", header, "GPL1261.matrix.",
    Sys.Date(), ".RData",
    sep = ""
    )
  )

#####
# Test
single.chip.enrichment.preRanked <- function (exprs, geneset)
# single chip enrichment based on median rank
# pre-ranked data 
{
    Ns <- ncol(exprs)
    gene.names <- rownames(exprs)
    geneset.names <- names(geneset)
    score.matrix <- matrix(0, nrow = length(geneset), ncol = Ns)
    for (i in 1:length(geneset)) {
        overlap <- intersect(geneset[[i]], gene.names)
        if (length(overlap) == 0) {
            score.matrix[i, ] <- NA
        }
        else {
              score.matrix[i, ] <- apply(exprs, 2, function(x) {
                median(x[overlap])
                })
    		}
    	}
    colnames(score.matrix) <- colnames(exprs)
    rownames(score.matrix) <- geneset.names
    return(score.matrix)
}

load("~/database-fingerprint-scripts/unionGenesets/DPD.Union.Mm.gs")
# calculate 1st geneset in first 5 arrays
testEnrichment <- single.chip.enrichment.preRanked(t(SCE.matrix[1:5,]),DPD.Union.Mm.gs)
# compare to values calcualted directly
testEnrichmentDirect <- matrix(nrow = length(DPD.Union.Mm.gs), ncol = 5)
rownames(testEnrichmentDirect) <- names(DPD.Union.Mm.gs)
colnames(testEnrichmentDirect) <- rownames(SCE.matrix)[1:5]
for (i in 1:5){
	tempfile = paste("~/Documents/Projects/Fingerprinting/data/Fingerprints/UnionDPD/UnionDPD_",
					 rownames(SCE.matrix)[i],
					 sep = "")
	load(tempfile)
	testEnrichmentDirect[,i]<-temp1[[1]]$SCG
	}

head(testEnrichment)	
head(testEnrichmentDirect)
cor(testEnrichment, testEnrichmentDirect,use = "complete.obs")
# looks like this has worked!! :)
# While 'here' calculate Jiantao's signature




#######

# remove na cols
#SCE.matrix <- SCE.matrix[,-c(1238, 18991)]
#SCE.matrix <- SCE.matrix[-c(1718,6508),]
# remove na rows
SCE.matrix <- SCE.matrix[rowSums(is.na(SCE.matrix)) == 0,]


# calculate correlation matrix
# create line by line
# multiply by 10000000 to maintain integer values

# need to split into mouse and human
cor.matrix <- matrix(nrow = ncol(SCE.matrix),
					 ncol = ncol(SCE.matrix))
colnames(cor.matrix)<-colnames(SCE.matrix)
rownames(cor.matrix)<-colnames(SCE.matrix)
print(object.size(cor.matrix), units = "auto")
# test two human matrices
GPL570a.cor.matrix <- cor.matrix
GPL570b.cor.matrix <- cor.matrix


GPL570.cor.matrix <- cor.matrix
GPL1261.cor.matrix <- cor.matrix




SCE.matrix.GPL570 <- SCE.matrix[
	rownames(SCE.matrix) %in% probe.data$GSM[probe.data$GPL == "GPL570"]
			,]


SCE.matrix.GPL1261 <- SCE.matrix[
	rownames(SCE.matrix) %in% probe.data$GSM[probe.data$GPL == "GPL1261"]
			,]


rm(cor.matrix)

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