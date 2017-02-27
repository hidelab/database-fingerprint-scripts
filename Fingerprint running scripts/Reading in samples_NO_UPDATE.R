# Creating updated list of GEO references
#
# Author: Gabriel Altschuler
# Status: Updated
# Timestamp: 20111128
# Script to read list of GEO samples in the current fingerprint version to update the fingerprint analysis with GEO files from exisitng fingerprint

GEOsamples.new<-GEO.metadata.matrix$GSM[GEO.metadata.matrix$GPL %in% names(chipframe)]

# scramble sample names to avoid bias in processing time asscociated with the larger chips
# apportion into 12 or 14
print(paste(length(GEOsamples.new), "to be fingerprinted", sep = " ")) 
GEOsamples.rand <- sample(GEOsamples.new)
len <- length(GEOsamples.rand) %/% 14

sampleset<-vector("list", 14)
for (i in 1:length(sampleset)){
	sampleset[[i]] <- GEOsamples.rand[(1+((i-1)*len)):(i*len)]
	}
sampleset[[14]] <- GEOsamples.rand[(1+(13*len)):length(GEOsamples.rand)]

