# estimating p-value distribution for random permutations of fingeprint vectors and lengths
# Author: Gabriel Altschuler
# Timestamp: 20010512
# Status: In progress
# Script to analyze the fingerprint distance to get empirical p-value


source("/Users/GabrielAltschuler/Dropbox/fingerprint/scripts/gabriel functions.R")
# load fingerprint matrix
load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Data/POE.matrix.0.001.2011-04-06")

# construct dummy fingerprint consensus
lengths<-50
reps<-5
sd<-vector("list", lengths)
mean<-vector("list", lengths)
for (j in 2:lengths){
	sd.temp<-vector("numeric", reps)
	mean.temp<-vector("numeric", reps)
	for (i in 1:reps){
		fingerprint<-rep(0,nrow(POE.matrix.0.001))
		n<-sample(0:j,1)
		values<-c(rep(-1,n), rep(1, (j - n)))
		fingerprint[sample(1:length(fingerprint), j)]<-values
		fingerprint.distance<-consensusDistance(as.data.frame(fingerprint), POE.matrix.0.001)
		sd.temp[i]<-sd(fingerprint.distance)
		mean.temp[i]<-mean(fingerprint.distance)
		}
	sd[[j]]<-sd.temp
	mean[[j]]<-mean.temp
	}
	
	
# now need to plot out

plot(1:49, lapply(mean[-1], mean))
plot(1:49, lapply(sd[-1], mean))

# looks like mean tends to 0.5
