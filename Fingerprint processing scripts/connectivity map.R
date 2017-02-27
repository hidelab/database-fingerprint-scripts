# Pathway connectivity map
# Author: Gabriel Altschuler
# Timestamp: 20110426
# Status: Complete
# Script to extract the fingerprints of the CMAP data

library(GEOmetadb)
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
gsm.char = dbGetQuery(con, "select gsm, characteristics_ch1 from gsm")
gsm.source = dbGetQuery(con, "select gsm, source_name_ch1 from gsm")
connectivityMap<-geoConvert("GSE5258", c("gsm"))$gsm
dbDisconnect(con)

# load experimental matrix
load("/data/shared/Fingerprint/Experimental.POE.matrix.2011-04-15")
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/Analysis")

# all cmap arrays have been fingerprinted
sum(!(connectivityMap[,2] %in% colnames(POE.matrix)))
Experimental.connectivityMap.POE<-POE.matrix[,connectivityMap[,2]]

# load full matrix
defineFile<-function(){
readline("enter filename: ")
}
dataPath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/"
dir(path = dataPath,
    pattern = "POE.matrix.20")
print("select POE file")
POEfile<-defineFile()
load(paste(dataPath, POEfile, sep = ""))

connectivityMap.POE<-POE.matrix[,connectivityMap[,2]]

# save to shared directory
save(Experimental.connectivityMap.POE, file = "/data/shared/Fingerprint/Experimental.connectivityMap.POE.RData")
save(connectivityMap.POE, file = "/data/shared/Fingerprint/connectivityMap.POE.RData")

# Threshold and save
connectivityMap.POE.001<-(connectivityMap.POE > 0.001) - (connectivityMap.POE < 0.001)
Experimental.connectivityMap.POE.001<-(Experimental.connectivityMap.POE > 0.001) - (Experimental.connectivityMap.POE < 0.001)
save(Experimental.connectivityMap.POE.001, file = "/data/shared/Fingerprint/Experimental.connectivityMap.POE.001.RData")
save(connectivityMap.POE.001, file = "/data/shared/Fingerprint/connectivityMap.POE.001.RData")

# Save metadata
connectivityMap.meta<-gsm.char[match(connectivityMap[,2], gsm.char[,1]),]
save(connectivityMap.meta, file = "/data/shared/Fingerprint/connectivityMap.meta.RData")

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/Analysis")
pdf("CMAP.POE.pdf")

# Just out of interest, plot histograms for the experimental datasets

signedlog2<-function(x){sign(x)*log((1+1000*abs(x)), 2)}

pdf("/home/galtschu2/Documents/Projects/Fingerprinting/Analysis/Experimental.CMAP.histograms.pdf")
for (i in 1:nrow(Experimental.connectivityMap.POE)){
  temp<-hist(signedlog2(Experimental.connectivityMap.POE[i,]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  plot(x = temp$mids, y = temp$counts, type = "h", lwd = 2.5, xlab= "Geneset POE score", xaxt = "n", ylab = "frequency", main = rownames(Experimental.connectivityMap.POE)[i])
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  }
dev.off()

# End of script