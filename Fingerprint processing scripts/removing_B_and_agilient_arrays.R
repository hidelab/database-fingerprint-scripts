# removing B arrays from the pathway fingerprint
# Author: Gabriel Altschuler
# Timestamp: 08022011
# short script to remove the B and agilent arrays from the pathprint dataframe and chipframe
library(pathprint.v0.3.beta4)
data(GEO.metadata.matrix)
data(GEO.fingerprint.matrix)
GEO.metadata.matrix<-GEO.metadata.matrix[
  !(GEO.metadata.matrix$GPL %in% c("GPL340", "GPL97")),]

GEO.fingerprint.matrix<-GEO.fingerprint.matrix[,
  colnames(GEO.fingerprint.matrix) %in% GEO.metadata.matrix$GSM]

GEO.metadata.matrix<-GEO.metadata.matrix[
  GEO.metadata.matrix$GSM %in% colnames(GEO.fingerprint.matrix),]

save(GEO.metadata.matrix,
     file = "/Users/GabrielAltschuler/Dropbox/fingerprint/data/GEO.metadata.matrix.RData")

save(GEO.fingerprint.matrix,
     file = "/Users/GabrielAltschuler/Dropbox/fingerprint/data/GEO.fingerprint.matrix.RData")

# also remove these arrays from the chipframe     
agilents<-c("GPL1708", "GPL4133", "GPL4134", "GPL6466", "GPL6480", "GPL6848", "GPL7202", "GPL887",  "GPL891")

badArrays<-c(agilents, "GPL340", "GPL97")

load("/Users/GabrielAltschuler/Dropbox/fingerprint/data/chipframe.RData")

chipframe<-chipframe[!(names(chipframe) %in% badArrays)]

save(chipframe, file = "/Users/GabrielAltschuler/Dropbox/fingerprint/data/chipframe.RData")