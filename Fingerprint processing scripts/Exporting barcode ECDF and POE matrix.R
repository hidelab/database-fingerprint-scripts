# Producing barcode matrix from SCG.frame and SCG.POE
# load data frame containing platform data and SCG.frame
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/platform_SCG_frames.RData")
# load SCG.POE
setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
load("data/SCG.POE.RData")

# now want to output the barcode matrix files

barcode.celltypes<-read.delim("Barcode/barcode_figure2_data.txt", sep = "\t", stringsAsFactors = FALSE)

# check that arrays are present in matrix

barcode.celltypes$DB_ID[!(barcode.celltypes$DB_ID %in% colnames(t(SCG.frame)))]

# have 541 of 552 samples in matrix - the missing ones are the yeast arrays and one stray lung sample
sum((barcode.celltypes$DB_ID %in% colnames(t(SCG.frame))))

barcode.SCG.matrix<-(t(SCG.frame))[,colnames(t(SCG.frame)) %in% barcode.celltypes$DB_ID]
barcode.POE.matrix<-t(SCG.POE)[,colnames(t(SCG.POE)) %in% barcode.celltypes$DB_ID]

## ERROR IN THIS PART - see below
# threshold at 5,10,15,20,25,30% for SCG at 
threshold<-c(5,10,15,20,25,30)
for (i in 1:length(threshold)){
  high<-apply(t(SCG.frame),1, quantile, probs = (100-threshold[i])/100) 
  low<-apply(t(SCG.frame),1, quantile, probs = threshold[i]/100)
	barcode.ECDF.matrix.threshold<-((barcode.SCG.matrix)>high)-((barcode.SCG.matrix)<low)
	barcode.ECDF.matrix.threshold[is.na(barcode.ECDF.matrix.threshold)]<-0
	assign(paste("barcode.ECDF.matrix.threshold.", threshold[i], sep = ""), barcode.ECDF.matrix.threshold)
	save(	list = paste("barcode.ECDF.matrix.threshold.", threshold[i], sep = ""),
			file = paste("Barcode/barcode.ECDF.matrix.threshold.", threshold[i], ".RData", sep = "")
			)
	print(threshold[i])
  }
  
# N.B. colnames are in a different order to the previous matrix
#output the barcode.POE.matrix as is - can be thresholded or used raw by Rebecca

###

########
# The ECDF normalization did not take into account the chip type!
# need to re-do this
platforms<-as.character(unique(platform.frame$Platform))
threshold<-c(5,10,15,20,25,30)
for (i in 1:length(threshold)){
  ECDF.matrix.threshold <- matrix(nrow = nrow(t(SCG.frame)), ncol = ncol(t(SCG.frame)))
  dimnames(ECDF.matrix.threshold) <- dimnames(t(SCG.frame))
  for (j in 1:length(platforms)){
		frame<-t(SCG.frame[grep(platforms[j], platform.frame$Platform), ])
		high<-apply(frame, 1, quantile, probs = (100-threshold[i])/100) 
		low<-apply(frame, 1, quantile, probs = threshold[i]/100)
		ECDF.matrix.threshold[, grep(platforms[j], platform.frame$Platform)]<-((frame)>high)-((frame)<low)
		ECDF.matrix.threshold[, grep(platforms[j], platform.frame$Platform)][is.na(ECDF.matrix.threshold[, grep(platforms[j], platform.frame$Platform)])]<-0
	}	
	barcode.ECDF.matrix.threshold<-ECDF.matrix.threshold[,colnames(ECDF.matrix.threshold) %in% barcode.celltypes$DB_ID]
	assign(paste("barcode.ECDF.matrix.threshold.", threshold[i], sep = ""), barcode.ECDF.matrix.threshold)
	save(	list = paste("barcode.ECDF.matrix.threshold.", threshold[i], sep = ""),
				file = paste("Barcode/barcode.ECDF.matrix.threshold.", threshold[i], ".RData", sep = "")
		)
  print(threshold[i])
  }





save(barcode.POE.matrix, file = "Barcode/barcode.POE.matrix.RData")


