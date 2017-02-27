# Optimization of threshold for pluripotency signature
# Procedure
# 1) Split ES/iPS set into 2
# 2) Construct consensus fingerprint based on set 1

# 3) Measure 2 parameters; a) median distance of full GEO set from ES/iPS signature b) median distance of ES/iPS set 2

# 4) Sample 100x
# 5) Repeat for various fingerprint thresholds
# Run on server

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data")
load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Analysis/optimizingThreshold.RData")
database<-read.delim("GEO_ES_iPS_non_stem.txt")
database.valid<-database[(database$GSM %in% colnames(ReNorm.matrix.ternary)),]
database.human<-database.valid[grep("Human", database.valid$Species),]
pluripotents<-as.character(database.human[c(grep("ES", database.human$SimpleCellType) ,grep("^iPS", database.human$SimpleCellType)),3])

##### Could run additional thresholding on server, already have 5, 10, 15, 20, 25
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix7_10_10.RData")

# re-apply thresholds, higher cutoff as within smaller groups
thresholdedMatrixSet<-vector("list", 0)
x<-1:20
for (i in 1:length(x)){
	high<-(100-x[i])/100
	low<-x[i]/100
	thresholdedMatrix<-((ReNorm.matrix)>high)-((ReNorm.matrix)<low)
	thresholdedMatrix[is.na(thresholdedMatrix)]<-0
	thresholdedMatrixSet<-append(thresholdedMatrixSet, list(thresholdedMatrix))
	}
names(thresholdedMatrixSet)<-paste("Threshold", x, sep = "_")
save(thresholdedMatrixSet, "/home/galtschu2/Documents/Projects/Fingerprinting/data/7_10_10_ThresoldedSet.RData")


setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data")
database<-read.delim("GEO_ES_iPS_non_stem.txt")
database.valid<-database[(database$GSM %in% colnames(thresholdedMatrixSet[[1]])),]
database.human<-database.valid[grep("Human", database.valid$Species),]
# remove outlier arrays
outliers<-c("GSM423940", "GSM423941", "GSM423942", "GSM423943", "GSM423944", "GSM423945")
database.human<-database.human[!(database.human$GSM %in% outliers),]

pluripotents<-as.character(database.human[c(grep("ES", database.human$SimpleCellType) ,grep("^iPS", database.human$SimpleCellType)),3])

pluripotentStat<-function(thresholded.fingerprints, signature.threshold){
	#signature.threshold<-0.75
	GEO.distance<-rep(0,100)
	GEO.sd<-rep(0,100)
	control.distance<-rep(0,100)
	control.sd<-rep(0,100)
	nPathways<-rep(0,100)
	Pathways<-vector("list", 100)
	for (i in 1:100){
		pluripotent.test.sample<-sample(pluripotents, round(0.5*length(pluripotents),0))
		pluripotent.control.sample<-pluripotents[!(pluripotents %in% pluripotent.test.sample)]

		pluripotent.test.fingerprints<-thresholded.fingerprints[,match(pluripotent.test.sample, colnames(thresholded.fingerprints))]
	
		pluripotent.test.signature<-rowMeans(pluripotent.test.fingerprints)
		pluripotent.test.signature.fingerprint<-rep(0,length(names(pluripotent.test.signature)))
		names(pluripotent.test.signature.fingerprint)<-names(pluripotent.test.signature)
		pluripotent.test.signature.fingerprint[pluripotent.test.signature>signature.threshold]<-1
		pluripotent.test.signature.fingerprint[pluripotent.test.signature<(-signature.threshold)]<-(-1)

		# build full distribution
		pluripotent.difference<-thresholded.fingerprints[(abs(pluripotent.test.signature)>signature.threshold),] - 		pluripotent.test.signature.fingerprint[(abs(pluripotent.test.signature)>signature.threshold)]
		nPathways[i]<-sum(abs(pluripotent.test.signature)>signature.threshold)
		pluripotent.distance<-colSums(abs(pluripotent.difference))/nPathways[i]
		GEO.distance[i]<-median(pluripotent.distance)
		GEO.sd[i]<-sd(pluripotent.distance)
		control.distance[i]<-median(pluripotent.distance[pluripotent.control.sample])
		control.sd[i]<-sd(pluripotent.distance[pluripotent.control.sample])
		Pathways[i]<-list(names(pluripotent.test.signature[(abs(pluripotent.test.signature)>signature.threshold)]))
		}
	return(list(GEO.distance = GEO.distance, GEO.sd = GEO.sd, control.distance = control.distance, control.sd = control.sd, nPathways = nPathways, Pathways = Pathways))
	}

# now expand to provide optimization over consensus threshold range as the thresholded fingerprints.
pluripotentStat.range<-function(thresholded.fingerprints, signature.threshold.range){
	#signature.threshold<-0.75
	GEO.distance<-matrix(ncol = length(signature.threshold.range), nrow = 100)
	GEO.sd<-matrix(ncol = length(signature.threshold.range), nrow = 100)
	control.distance<-matrix(ncol = length(signature.threshold.range), nrow = 100)
	control.sd<-matrix(ncol = length(signature.threshold.range), nrow = 100)
	nPathways<-matrix(ncol = length(signature.threshold.range), nrow = 100)
	colnames(GEO.distance)<-paste("Threshold", signature.threshold.range)
	colnames(GEO.sd)<-paste("Threshold", signature.threshold.range)
	colnames(control.distance)<-paste("Threshold", signature.threshold.range)
	colnames(control.sd)<-paste("Threshold", signature.threshold.range)
	colnames(nPathways)<-paste("Threshold", signature.threshold.range)
	
	#Pathways<-vector("list", 100)
	for (j in 1:length(signature.threshold.range)){
		for (i in 1:100){
			pluripotent.test.sample<-sample(pluripotents, round(0.5*length(pluripotents),0))
			pluripotent.control.sample<-pluripotents[!(pluripotents %in% pluripotent.test.sample)]

			pluripotent.test.fingerprints<-thresholded.fingerprints[,match(pluripotent.test.sample, colnames(thresholded.fingerprints))]
	
			pluripotent.test.signature<-rowMeans(pluripotent.test.fingerprints)
			pluripotent.test.signature.fingerprint<-rep(0,length(names(pluripotent.test.signature)))
			names(pluripotent.test.signature.fingerprint)<-names(pluripotent.test.signature)
			pluripotent.test.signature.fingerprint[pluripotent.test.signature>signature.threshold.range[j]]<-1
			pluripotent.test.signature.fingerprint[pluripotent.test.signature<(-signature.threshold.range[j])]<-(-1)

			# build full distribution
			pluripotent.difference<-thresholded.fingerprints[(abs(pluripotent.test.signature)>signature.threshold.range[j]),] - 			pluripotent.test.signature.fingerprint[(abs(pluripotent.test.signature)>signature.threshold.range[j])]
			nPathways[i,j]<-sum(abs(pluripotent.test.signature)>signature.threshold.range[j])
			pluripotent.distance<-colSums(abs(pluripotent.difference))/nPathways[i,j]
			GEO.distance[i,j]<-median(pluripotent.distance)
			GEO.sd[i,j]<-sd(pluripotent.distance)
			control.distance[i,j]<-median(pluripotent.distance[pluripotent.control.sample])
			control.sd[i,j]<-sd(pluripotent.distance[pluripotent.control.sample])
			#Pathways[i,j]<-list(names(pluripotent.test.signature[(abs(pluripotent.test.signature)>signature.threshold.range[j])]))
			}
		}
		
	return(list(GEO.distance = GEO.distance, GEO.sd = GEO.sd, control.distance = control.distance, control.sd = control.sd, nPathways = nPathways))
	}






x<-1:20
Threshold.review<-vector("list", length(x))
for (i in 1:length(x)){
	Threshold.review[i]<-list(pluripotentStat(thresholdedMatrixSet[[i]], signature.threshold = 0.75))
	}

save(Threshold.review, file = "Threshold.review.RData")

# expand to wider threshold - start with empty workspace so have to save and re-load bits as going along
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix7_10_10.RData")
thresholdedMatrixSet<-vector("list", 0)
x<-21:35
for (i in 1:length(x)){
	high<-(100-x[i])/100
	low<-x[i]/100
	thresholdedMatrix<-((ReNorm.matrix)>high)-((ReNorm.matrix)<low)
	thresholdedMatrix[is.na(thresholdedMatrix)]<-0
	thresholdedMatrixSet<-append(thresholdedMatrixSet, list(thresholdedMatrix))
	}
names(thresholdedMatrixSet)<-paste("Threshold", x, sep = "_")
save(thresholdedMatrixSet, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/7_10_10_ThresoldedSet_21_35.RData")


Threshold.review<-vector("list", length(x))
for (i in 1:length(x)){
	Threshold.review[i]<-list(pluripotentStat(thresholdedMatrixSet[[i]], signature.threshold = 0.75))
	}

save(Threshold.review, file = "Threshold.review21_35.RData")

###############
# repeat using the 2D matrix version
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix7_10_10.RData")

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data")
database<-read.delim("GEO_ES_iPS_non_stem.txt")
database.valid<-database[(database$GSM %in% colnames(ReNorm.matrix)),]
database.human<-database.valid[grep("Human", database.valid$Species),]
# remove outlier arrays
outliers<-c("GSM423940", "GSM423941", "GSM423942", "GSM423943", "GSM423944", "GSM423945")
database.human<-database.human[!(database.human$GSM %in% outliers),]

pluripotents<-as.character(database.human[c(grep("ES", database.human$SimpleCellType) ,grep("^iPS", database.human$SimpleCellType)),3])




x<-1:35
pb <- txtProgressBar(min = 0, max = length(x), style = 3)
Threshold.review<-vector("list", length(x))
for (i in 1:length(x)){
	high<-(100-x[i])/100
	low<-x[i]/100
	thresholdedMatrix<-((ReNorm.matrix)>high)-((ReNorm.matrix)<low)
	thresholdedMatrix[is.na(thresholdedMatrix)]<-0
	Threshold.review[i]<-list(pluripotentStat.range(thresholdedMatrix, signature.threshold.range = seq(0.5, 1, 0.05)))
	setTxtProgressBar(pb, i)
	}
save(Threshold.review, file = "Threshold.review_2D.RData")


# now on local

load("/Users/GabrielAltschuler/Documents/Projects/fingerprinting/data/Threshold.review.RData")
# need to add second part - 21-35
Threshold.review1<-Threshold.review
load("/Users/GabrielAltschuler/Documents/Projects/fingerprinting/data/Threshold.review21_35.RData")
Threshold.review2<-append(Threshold.review1, Threshold.review)
x<-1:35
pathwayNumber<-vector("list", length(x))
for (i in x){	
	pathwayNumber<-c(pathwayNumber, Threshold.review2[[i]]$nPathways)
	}
	
plot(unlist(pathwayNumber))

GEO.distance<-vector("list", length(x))
for (i in x){	
	GEO.distance<-c(GEO.distance, Threshold.review2[[i]]$GEO.distance)
	}

plot(unlist(GEO.distance))

control.distance<-vector("list", length(x))
for (i in x){	
	control.distance<-c(control.distance, Threshold.review2[[i]]$control.distance)
	}
nPathways<-vector("list", length(x))
for (i in x){	
	nPathways<-c(nPathways,Threshold.review2[[i]]$nPathways)
	}
GEO.sd<-vector("list", length(x))
for (i in x){	
	GEO.sd<-c(GEO.sd,Threshold.review2[[i]]$GEO.sd)
	}
control.sd<-vector("list", length(x))
for (i in x){	
	control.sd<-c(control.sd,Threshold.review2[[i]]$control.sd)
	}

plot(unlist(control.distance))

plot(unlist(GEO.distance)-unlist(control.distance))

names<-c(
		rep(1, 100),
		rep(2, 100),
		rep(3, 100),
		rep(4, 100),
		rep(5, 100),
		rep(6, 100),
		rep(7, 100),
		rep(8, 100),
		rep(9, 100),
		rep(10, 100),
		rep(11, 100),
		rep(12, 100),
		rep(13, 100),
		rep(14, 100),
		rep(15, 100),
		rep(16, 100),
		rep(17, 100),
		rep(18, 100),
		rep(19, 100),
		rep(20, 100),
		rep(21, 100),
		rep(22, 100),
		rep(23, 100),
		rep(24, 100),
		rep(25, 100),
		rep(26, 100),
		rep(27, 100),
		rep(28, 100),
		rep(29, 100),
		rep(30, 100),
		rep(31, 100),
		rep(32, 100),
		rep(33, 100),
		rep(34, 100),
		rep(35, 100)
		)
		
info<-data.frame(name = names, GEO.distance = unlist(GEO.distance), control.distance = unlist(control.distance), nPathways = unlist(nPathways), GEO.sd = unlist(GEO.sd), control.sd = unlist(control.sd))

info[,7]<-((info$GEO.distance-info$GEO.sd)-(info$control.distance+info$control.sd))
info[,8]<-((info$GEO.distance-2*info$GEO.sd)-(info$control.distance+2*info$control.sd))
info[,9]<-((info$GEO.distance-3*info$GEO.sd)-(info$control.distance+3*info$control.sd))
setwd("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Analysis")

pdf("Threshold summary.pdf")
boxplot(info$GEO.distance ~ info$name, xlab = "Threshold", ylab = "Normalized distance", main = "Corpus median")
boxplot(info$control.distance ~ info$name, xlab = "Threshold", ylab = "Normalized distance", main = "ES/iPS test median")
boxplot(info$nPathways ~ info$name, xlab = "Threshold", ylab = "Number of pathways", main = "Composition of pluripotency fingerprint")
boxplot(info$GEO.sd ~ info$name, xlab = "Threshold", ylab = "Normalized distance", main = "Corpus stdev")
boxplot(info$control.sd ~ info$name, xlab = "Threshold", ylab = "Normalized distance", main = "ES/iPS test stdev")
boxplot(info[,7] ~ info$name, xlab = "Threshold", ylab = "Normalized distance", main = "Spread")
boxplot(info[,8] ~ info$name, xlab = "Threshold", ylab = "Normalized distance", main = "Spread - 2 SD")
boxplot(info[,9] ~ info$name, xlab = "Threshold", ylab = "Normalized distance", main = "Spread - 3 SD")
dev.off()




save.image("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Analysis/optimizingThreshold.RData")



###### 2D - optimization
load("/Users/GabrielAltschuler/Documents/Projects/fingerprinting/data/Threshold.review_2D.RData")
x<-1:35
pathwayNumberMatrix<-matrix(ncol = 11, nrow = length(x))
for (i in x){	
	pathwayNumberMatrix[i,]<-apply(Threshold.review[i][[1]]$nPathways, 2, median)
	}


GEODistanceMatrix<-matrix(ncol = 11, nrow = length(x))
for (i in x){	
	GEODistanceMatrix[i,]<-apply(Threshold.review[i][[1]]$GEO.distance, 2, median)
	}

ControlDistanceMatrix<-matrix(ncol = 11, nrow = length(x))
for (i in x){	
	ControlDistanceMatrix[i,]<-apply(Threshold.review[i][[1]]$control.distance, 2, median)
	}



GEOSDMatrix<-matrix(ncol = 11, nrow = length(x))
for (i in x){	
	GEOSDMatrix[i,]<-apply(Threshold.review[i][[1]]$GEO.sd, 2, median)
	}



ControlSDMatrix<-matrix(ncol = 11, nrow = length(x))
for (i in x){	
	ControlSDMatrix[i,]<-apply(Threshold.review[i][[1]]$control.sd, 2, median)
	}


colors <- colorRampPalette(c("green", "black", "black", "red"))
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"), bias = 0.5)


image(1:35, seq(0.5, 0.95, 0.05), pathwayNumberMatrix[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = jet.colors(100), main = "nPathways")
contour(1:35, seq(0.5, 0.95, 0.05), pathwayNumberMatrix[,1:10], add = TRUE, col = "peru")
image(1:35, seq(0.5, 0.95, 0.05), GEODistanceMatrix[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = colors(100), main = "GEO distance")
contour(1:35, seq(0.5, 0.95, 0.05), GEODistanceMatrix[,1:10], add = TRUE, col = "peru")
image(1:35, seq(0.5, 0.95, 0.05), ControlDistanceMatrix[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = colors(100), main = "Control distance")
contour(1:35, seq(0.5, 0.95, 0.05), ControlDistanceMatrix[,1:10], add = TRUE, col = "peru")
image(1:35, seq(0.5, 0.95, 0.05), GEOSDMatrix[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = colors(100), main = "GEO SD")
contour(1:35, seq(0.5, 0.95, 0.05), GEOSDMatrix[,1:10], add = TRUE, col = "peru")
image(1:35, seq(0.5, 0.95, 0.05), ControlSDMatrix[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = colors(100), main = "Control SD")
contour(1:35, seq(0.5, 0.95, 0.05), ControlSDMatrix[,1:10], add = TRUE, col = "peru")

diff.distance.1<-(GEODistanceMatrix-GEOSDMatrix)-(ControlDistanceMatrix+ControlSDMatrix)
diff.distance.2<-(GEODistanceMatrix-2*GEOSDMatrix)-(ControlDistanceMatrix+2*ControlSDMatrix)
diff.distance.3<-(GEODistanceMatrix-3*GEOSDMatrix)-(ControlDistanceMatrix+3*ControlSDMatrix)

image(1:35, seq(0.5, 0.95, 0.05), diff.distance.1[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = jet.colors(100), main = "diff.distance.1")
contour(1:35, seq(0.5, 0.95, 0.05), diff.distance.1[,1:10], add = TRUE, col = "peru", nlevels = 20)
image(1:35, seq(0.5, 0.95, 0.05), diff.distance.2[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = jet.colors(100), main = "diff.distance.2")
contour(1:35, seq(0.5, 0.95, 0.05), diff.distance.2[,1:10], add = TRUE, col = "peru", nlevels = 40)
image(1:35, seq(0.5, 0.95, 0.05), diff.distance.3[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = jet.colors(100), main = "diff.distance.3")
contour(1:35, seq(0.5, 0.95, 0.05), diff.distance.3[,1:10], add = TRUE, col = "peru", nlevels = 20)



max.dist<-max(diff.distance.3[,1:10])

maxConsensusThreshold<-(grep(max(diff.distance.3[,1:10]), diff.distance.3[,1:10]) %/% length(x)+1)
maxFingerprintThreshold<-grep(max(diff.distance.3[,1:10]), diff.distance.3[,1:10]) %% length(x)

maxFingerprintThreshold.data<-((Threshold.review[maxFingerprintThreshold][[1]]$GEO.distance-3*Threshold.review[maxFingerprintThreshold][[1]]$GEO.sd)-(Threshold.review[maxFingerprintThreshold][[1]]$control.distance+3*Threshold.review[maxFingerprintThreshold][[1]]$control.sd))[,1:10]



maxConsensusThreshold.data<-matrix(nrow = 100, ncol = 35)
colnames(maxConsensusThreshold.data)<-paste(11:45, "%", sep = "_")
for (i in 1:ncol(maxConsensusThreshold.data)){
	maxConsensusThreshold.data[,i]<-((Threshold.review[i][[1]]$GEO.distance[, maxConsensusThreshold]-3*Threshold.review[i][[1]]$GEO.sd[, maxConsensusThreshold])-(Threshold.review[i][[1]]$control.distance[, maxConsensusThreshold]+3*Threshold.review[i][[1]]$control.sd[, maxConsensusThreshold]))
		}



pdf("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Analysis/Threshold_2D_summary.pdf")
par(mfrow = c(2,2), lwd = 0.5)
par(cex = 0.5)

image(1:35, seq (50,95,5), diff.distance.3[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold (%)", col = jet.colors(100), main = "ES/iPS vs GEO corpus distance parameter")
contour(1:35, seq (50,95,5), diff.distance.3[,1:10], add = TRUE, col = "black", nlevels = 10)
abline (v = maxFingerprintThreshold, lty = 2)
abline (h = 75, lty = 2)

image(1:35, seq (50,95,5), pathwayNumberMatrix[,1:10], xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold (%)", col = jet.colors(100), main = "Number of pathways in ES/iPS fingerprint", cex = 0.5)
contour(1:35, seq (50,95,5), pathwayNumberMatrix[,1:10], add = TRUE, col = "black", nlevels = 10)
abline (v = 15, lty = 2)
abline (h = 75, lty = 2)



boxplot(t(maxFingerprintThreshold.data)~rownames(t(maxFingerprintThreshold.data)), las = 2, main = paste("Fingerprint threshold =", maxFingerprintThreshold, "%", sep = " "), names = seq (50,95,5), xlab = "Consensus threshold (%)", ylab = "Distance parameter")
boxplot(t(maxConsensusThreshold.data)~rownames(t(maxConsensusThreshold.data)), las = 2, levels = 1:35, main = paste("Consensus threshold =", seq(50,95,5)[maxConsensusThreshold], "%", sep = " "), names = 1:35, xlab = "Fingerprint Threshold (%)", ylab = "Distance parameter")
save(Threshold.review, file = "Threshold.review_2D.RData")

#### re-do 2D optimization on server using KS-test and new scripts

x<-1:35 # range for fingerprint threshold
y<-seq(0.5, 0.95, 0.05) # range for consensus threshold
z<-1:100 # number of times to sample data
pb <- txtProgressBar(min = 0, max = length(x), style = 3)
KS.review.pvalue<-array(dim = c(length(x), length(y), length(z)))
KS.review.statistic<-array(dim = c(length(x), length(y), length(z)))
KS.review.nPathways<-array(dim = c(length(x), length(y), length(z)))
for (i in 1:length(x)){
	high<-(100-x[i])/100
	low<-x[i]/100
	thresholdedMatrix<-((ReNorm.matrix)>high)-((ReNorm.matrix)<low)
	thresholdedMatrix[is.na(thresholdedMatrix)]<-0
		for (k in 1:length(z)){
			for (j in 1:length(y)){
				temp<-sampleKSdist(fingerprintframe = thresholdedMatrix,
											threshold = y[j],
											samplenames = pluripotents,
											proportion = 0.5
											)
				KS.review.pvalue[i,j,k]<-temp$pvalue
				KS.review.statistic[i,j,k]<-temp$statistic
				KS.review.nPathways[i,j,k]<-temp$nPathways
				}
			}
	setTxtProgressBar(pb, i)
	}
	
KS.Threshold.review = list(KS.review.pvalue, KS.review.statistic, KS.review.nPathways)
save(KS.Threshold.review, file = "KS.Threshold.review.RData")


######## back on local
load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/data/KS.Threshold.review.RData")
names(KS.Threshold.review)<-c("pvalue","statistic", "nPathways") 

# Plot image for the number of pathways in the fingerprint
image(1:35, seq(0.5, 0.95, 0.05), KS.Threshold.review$nPathways, xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = jet.colors(100), main = "nPathways")

# Plot image for the KS statistic that separates the distributions
image(1:35, seq(0.5, 0.95, 0.05), KS.Threshold.review$statistic, xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = jet.colors(100), main = "Statistic")

# Plot image for the KS pvalue that separates the distributions
image(1:35, seq(0.5, 0.95, 0.05), -log(KS.Threshold.review$pvalue,10), xlab = "Fingerprint Threshold (%)", ylab = "Median consensus Threshold", col = jet.colors(100), main = "-log(p-value)")

save.image("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Analysis/Threshold.2D.RData")
dev.off()