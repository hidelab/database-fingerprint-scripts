# Analysis of new experimental data in fingerprint system
# Timestamp: 20110419
# Author: Gabriel Altschuler
# Status: In progress
# Strategy
# 1) Normalize fingerprint as before
# 2) Sample arrays at stringent threshold level (5%)
# 3) Correlate with other pathways
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/Analysis")
load("Experimental.POE.Analysis.RData")

# load POE experimental matrix

load("/data/shared/Fingerprint/Experimental.POE.matrix.2011-04-15")


# plot histograms of experimental data

pdf("Experimental_POE_histograms.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(POE.matrix[i,], xlim = c(-1,1), breaks = seq(-1, 1, 0.001), xlab = "Geneset POE", main = rownames(POE.matrix)[i])
  plot(x = temp$mids, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, xlim = c(-1,1), ylab = "1+log(frequency)") 
  }
dev.off()

# Which arrays show high Langmoen up and low Langmoen down? 
# "LangmoenGBM_Up" "LangmoenGBM_Down" 

colnames(POE.matrix)[intersect(which(POE.matrix["LangmoenGBM_Up",]>0.95), which(POE.matrix["LangmoenGBM_Down",]<0.95))]

# which arrays show high nyeso geneset expresssion?
colnames(POE.matrix)[which(POE.matrix["nyeso1_ip_shortlist_ez",]>0.95)]

# Survey different cell types

CellMontageAffyArrays<-read.delim("/data/shared/Fingerprint/CellMontageAffyArrays.txt", stringsAsFactors = FALSE)
pluripotents<-CellMontageAffyArrays$GSM[CellMontageAffyArrays$Source == "PluripotentArrays"]
normal<-CellMontageAffyArrays$GSM[CellMontageAffyArrays$Type == "Normal"]
pluri_oligo<-read.delim("/data/shared/Fingerprint/GEO_ES_iPS_non_stem.txt", stringsAsFactors = FALSE)
outliers<-c("GSM423940", "GSM423941", "GSM423942", "GSM423943", "GSM423944", "GSM423945")
pluri_oligo<-pluri_oligo[!(pluri_oligo$GSM %in% outliers),]
oligopotent<-pluri_oligo$GSM[grep("Non-stem", pluri_oligo$SimpleCellType)]


library(GEOmetadb)
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
gsm.char = dbGetQuery(con, "select gsm, characteristics_ch1 from gsm")
gsm.source = dbGetQuery(con, "select gsm, source_name_ch1 from gsm")
dbDisconnect(con)

# Warnings are due to some issue with certain characteristics strings
# e.g. grep("test", gsm.char[97896,2])

HeLa<-gsm.char[grep("HeLa", gsm.char[,2], ignore.case = TRUE),1]
Jurkat<-gsm.char[grep("Jurkat", gsm.char[,2], ignore.case = TRUE),1]
MCF.7<-gsm.char[grep("MCF-7", gsm.char[,2], ignore.case = TRUE),1]
HEK.293<-gsm.char[grep("HEK293", gsm.char[,2], ignore.case = TRUE),1]
iPS<-gsm.char[grep("iPS", gsm.char[,2], ignore.case = TRUE),1]
U2OS<-gsm.char[grep("U2OS", gsm.char[,2], ignore.case = TRUE),1]
Caco.2<-gsm.char[grep("Caco-2", gsm.char[,2]),1]
HT29<-gsm.char[grep("HT29", gsm.char[,2], ignore.case = TRUE),1]
T.47D<-gsm.char[grep("T47D", gsm.char[,2], ignore.case = TRUE),1]
Mm.3T3<-gsm.char[grep("3T3", gsm.char[,2], ignore.case = TRUE),1]
Melanoma<-gsm.char[grep("Melanoma", gsm.char[,2], ignore.case = TRUE),1]
Fibroblast<-union(gsm.char[grep("^Fibroblast$", gsm.char[,2], ignore.case = TRUE),1],
                  gsm.char[grep("^Fibroblasts$", gsm.char[,2], ignore.case = TRUE),1]
                  )



Type = c(
rep("HeLa", sum(HeLa %in% colnames(POE.matrix))),
rep("Jurkat", sum(Jurkat %in% colnames(POE.matrix))),
rep("MCF.7", sum(MCF.7 %in% colnames(POE.matrix))),
rep("HEK.293", sum(HEK.293 %in% colnames(POE.matrix))),
#rep("iPS", sum(iPS %in% colnames(POE.matrix))),
rep("U2OS", sum(U2OS %in% colnames(POE.matrix))),
rep("Caco.2", sum(Caco.2 %in% colnames(POE.matrix))),
rep("HT29", sum(HT29 %in% colnames(POE.matrix))),
rep("ES/iPS", sum(pluripotents %in% colnames(POE.matrix))),
#rep("3T3", sum(Mm.3T3 %in% colnames(POE.matrix))),
rep("T-47D", sum(T.47D %in% colnames(POE.matrix))),
rep("Melanoma", sum(Melanoma %in% colnames(POE.matrix))),
rep("Fibroblast", sum(Fibroblast %in% colnames(POE.matrix)))
)

frame.GBM<-as.data.frame(matrix(ncol = 1 + nrow(POE.matrix), nrow = length(Type)))
colnames(frame.GBM)<- c("Type", rownames(POE.matrix))
frame.GBM$Type<-Type

for (i in 1:nrow(POE.matrix)){
  frame.GBM[,(i+1)]<-c(
  POE.matrix[i, HeLa[HeLa %in% colnames(POE.matrix)]],
  POE.matrix[i, Jurkat[Jurkat %in% colnames(POE.matrix)]],
  POE.matrix[i, MCF.7[MCF.7 %in% colnames(POE.matrix)]],
  POE.matrix[i, HEK.293[HEK.293 %in% colnames(POE.matrix)]],
  #POE.matrix[i, iPS[iPS %in% colnames(POE.matrix)]],
  POE.matrix[i, U2OS[U2OS %in% colnames(POE.matrix)]],
  POE.matrix[i, Caco.2[Caco.2 %in% colnames(POE.matrix)]],
  POE.matrix[i, HT29[HT29 %in% colnames(POE.matrix)]],
  POE.matrix[i, pluripotents[pluripotents %in% colnames(POE.matrix)]],
  #POE.matrix[i, Mm.3T3[Mm.3T3 %in% colnames(POE.matrix)]],
  POE.matrix[i, T.47D[T.47D %in% colnames(POE.matrix)]],
  POE.matrix[i, Melanoma[Melanoma %in% colnames(POE.matrix)]],
  POE.matrix[i, Fibroblast[Fibroblast %in% colnames(POE.matrix)]]
  )
  }
  
# plot
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/Analysis")
library(beanplot)
pdf("experimental_genesets_beanplots.pdf")
for (i in 1:nrow(POE.matrix)){
  beanplot(frame.GBM[,i+1]~frame.GBM[,1], las = 2, main = rownames(POE.matrix)[i])
  }
dev.off()

pdf("experimental_genesets_beanplots_lim.pdf")
for (i in 1:nrow(POE.matrix)){
  beanplot(frame.GBM[,i+1]~frame.GBM[,1], las = 2, main = rownames(POE.matrix)[i], ylim = c(-0.1, 0.1))
  }
dev.off()

pdf("experimental_genesets_beanplots_log.pdf")
for (i in 1:nrow(POE.matrix)){
  beanplot(log(200+frame.GBM[,i+1],2)~frame.GBM[,1], las = 2, main = rownames(POE.matrix)[i], overallline = "median")
  }
dev.off()

# Try plotting the melanoma samples CD133 hi vs low
pdf("CD133hivsCD133low.pdf")
plot(frame.GBM$CD133_common_Neg_peaks ~ frame.GBM$CD133_common_Pos_peaks, ylab = "CD133_common_Neg_peaks", xlab = "CD133_common_Pos_peaks", xlim = c(-1,1), ylim = c(-1,1))
plot(frame.GBM$CD133_common_Neg_peaks[frame.GBM$Type == "Melanoma"] ~ frame.GBM$CD133_common_Pos_peaks[frame.GBM$Type == "Melanoma"],
    main = "Melanoma samples only", ylab = "CD133_common_Neg_peaks", xlab = "CD133_common_Pos", xlim = c(-1,1), ylim = c(-1,1))
plot(frame.GBM$cd133hi_downreg ~ frame.GBM$cd133hi_upreg, ylab = "cd133hi_downreg", xlab = "cd133hi_upreg", xlim = c(-1,1), ylim = c(-1,1))
plot(frame.GBM$cd133hi_downreg[frame.GBM$Type == "Melanoma"] ~ frame.GBM$cd133hi_upreg[frame.GBM$Type == "Melanoma"],
    main = "Melanoma samples only", ylab = "cd133hi_downreg", xlab = "cd133hi_upreg", xlim = c(-1,1), ylim = c(-1,1))

dev.off()


pdf("Experimental_POE_histograms_with_pluripotents_twosided.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(POE.matrix[i,], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.pluripotents<-hist(POE.matrix[i,pluripotents], breaks = seq(-1, 1, 0.001), plot = FALSE)
  par(mfcol=c(1,2))
  # plot negative values
  plot(x = -temp$mids, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlim = c(1, 0.0005), xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)", log = "x") 
  points(x = -temp.pluripotents$mids, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "red", lwd = 2.5)
  # plot positive values
  plot(x = temp$mids, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlim = c(0.0005, 1), xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)", log = "x") 
  points(x = temp.pluripotents$mids, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "red", lwd = 2.5)
  }
dev.off()

pdf("Experimental_POE_histograms_with_pluripotents.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(POE.matrix[i,], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.pluripotents<-hist(POE.matrix[i,pluripotents], breaks = seq(-1, 1, 0.001), plot = FALSE)
  xvals<-temp$mids
  xvals[temp$mids>0]<-log(temp$mids[temp$mids>0],2)-log(0.00025,2)
  xvals[temp$mids<0]<-(-log(-temp$mids[temp$mids<0],2)+log(0.00025,2))
  # this creates two peaks at zero, need to offset centers
  #xvals[1000]<-(-log)
  #xvals[1001]<-0.1
  plot(x = xvals, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)") 
  points(x = xvals, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "red", lwd = 2.5)
  }
dev.off()

pdf("Experimental_POE_histograms_with_Melanoma.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(POE.matrix[i,], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.Melanoma<-hist(POE.matrix[i,Melanoma[Melanoma %in% colnames(POE.matrix)]], breaks = seq(-1, 1, 0.001), plot = FALSE)
  xvals<-temp$mids
  xvals[temp$mids>0]<-log(temp$mids[temp$mids>0],2)-log(0.00025,2)
  xvals[temp$mids<0]<-(-log(-temp$mids[temp$mids<0],2)+log(0.00025,2))
  # this creates two peaks at zero, need to offset centers
  #xvals[1000]<-(-log)
  #xvals[1001]<-0.1
  plot(x = xvals, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)") 
  points(x = xvals, y = 1+log(temp.Melanoma$counts,2), type = "h", col = "red", lwd = 2.5)
  }
dev.off()

pdf("Experimental_POE_histograms_with_pluripotents_and_fibroblasts.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(POE.matrix[i,], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.pluripotents<-hist(POE.matrix[i,pluripotents], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.fibroblasts<-hist(POE.matrix[i,Fibroblast[Fibroblast %in% colnames(POE.matrix)]], breaks = seq(-1, 1, 0.001), plot = FALSE)
  xvals<-temp$mids
  xvals[temp$mids>0]<-log(temp$mids[temp$mids>0],2)-log(0.00025,2)
  xvals[temp$mids<0]<-(-log(-temp$mids[temp$mids<0],2)+log(0.00025,2))
  # this creates two peaks at zero, need to offset centers  #xvals[1000]<-(-log)

  #xvals[1001]<-0.1
  plot(x = xvals, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)") 
  points(x = xvals, y = 1+log(temp.fibroblasts$counts,2), type = "h", col = "blue", lwd = 2.5)
  points(x = xvals, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "red", lwd = 2.5)
  }
dev.off()

pdf("Experimental_POE_histograms_with_pluripotents_and_normal.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(POE.matrix[i,], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.pluripotents<-hist(POE.matrix[i,pluripotents], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.normal<-hist(POE.matrix[i,normal[normal %in% colnames(POE.matrix)]], breaks = seq(-1, 1, 0.001), plot = FALSE)
  xvals<-temp$mids
  xvals[temp$mids>0]<-log(temp$mids[temp$mids>0],2)-log(0.00025,2)
  xvals[temp$mids<0]<-(-log(-temp$mids[temp$mids<0],2)+log(0.00025,2))
  # this creates two peaks at zero, need to offset centers
  #xvals[1000]<-(-log)
  #xvals[1001]<-0.1
  plot(x = xvals, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)") 
  points(x = xvals, y = 1+log(temp.normal$counts,2), type = "h", col = "blue", lwd = 2.5)
  points(x = xvals, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "red", lwd = 2.5)
  }
dev.off()

pdf("Experimental_POE_histograms_with_pluripotents_and_oligopotent.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(POE.matrix[i,], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.pluripotents<-hist(POE.matrix[i,pluripotents], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.oligopotent<-hist(POE.matrix[i,oligopotent[oligopotent %in% colnames(POE.matrix)]], breaks = seq(-1, 1, 0.001), plot = FALSE)
  xvals<-temp$mids
  xvals[temp$mids>0]<-log(temp$mids[temp$mids>0],2)-log(0.00025,2)
  xvals[temp$mids<0]<-(-log(-temp$mids[temp$mids<0],2)+log(0.00025,2))
  # this creates two peaks at zero, need to offset centers
  #xvals[1000]<-(-log)
  #xvals[1001]<-0.1
  plot(x = xvals, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)") 
  points(x = xvals, y = 1+log(temp.oligopotent$counts,2), type = "h", col = "blue", lwd = 2.5)
  points(x = xvals, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "#FF000050", lwd = 2.5)
  }
dev.off()

# Stop here - to do items include adding more fibroblast samples and possibly curating up one of the other cell lines (non-stemmy cancer cell line..?)

# check the list of fibroblast cells and add the oligopotents from the stem cell set

Fibroblasts.curated<-c("GSM248201", "GSM248202", "GSM248204", "GSM248209", "GSM248210", "GSM248213", "GSM248214", "GSM249027")
Fibroblasts.curated<-union(Fibroblasts.curated, oligopotent)
Fibroblasts.curated<-Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]

# curating additional fibroblasts
Fibroblast.additional<-union( gsm.source[grep("^Fibroblast$", gsm.source[,2], ignore.case = TRUE),1],
                              gsm.source[grep("^Fibroblasts$", gsm.source[,2], ignore.case = TRUE),1]
                              )

Fibroblast.additional<-Fibroblast.additional[Fibroblast.additional %in% colnames(POE.matrix)]

# found new ips experiment, GSE18111/GSE18226, not yet curated into iPS set, and perhaps GSE2375
"GSM452727" %in% pluripotents

Fibroblast.additional<-c(Fibroblast.additional, c("GSM467903","GSM467907","GSM467908"))
Fibroblast.additional<-setdiff(Fibroblast.additional, c("GSM629550","GSM629551", "GSM629552", "GSM628667", "GSM628668", "GSM628669", "GSM628670", "GSM628671", "GSM628672","GSM628673", "GSM628674", "GSM628675", "GSM628676"))

Fibroblasts.curated<-union(Fibroblasts.curated, Fibroblast.additional)

Fibroblast.additional.2<-union( gsm.source[grep("Fibroblast", gsm.source[,2], ignore.case = TRUE),1],
                              gsm.source[grep("Fibroblasts", gsm.source[,2], ignore.case = TRUE),1]
                              )
Fibroblast.additional.2<-Fibroblast.additional.2[Fibroblast.additional.2 %in% colnames(POE.matrix)]
Fibroblast.additional.2<-setdiff(Fibroblast.additional.2, union( gsm.source[grep("^Fibroblast$", gsm.source[,2], ignore.case = TRUE),1],gsm.source[grep("^Fibroblasts$", gsm.source[,2], ignore.case = TRUE),1]))

#### curated upto 300

Fibroblasts.curated.2<-readLines("/data/shared/Fingerprint/curatedCellTypes/Fibroblasts.txt")
Fibroblasts.curated<-union(Fibroblasts.curated, Fibroblasts.curated.2)

pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(POE.matrix[i,], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.pluripotents<-hist(POE.matrix[i,pluripotents], breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.Fibroblasts.curated<-hist(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]], breaks = seq(-1, 1, 0.001), plot = FALSE)
  xvals<-temp$mids
  xvals[temp$mids>0]<-log(temp$mids[temp$mids>0],2)-log(0.00025,2)
  xvals[temp$mids<0]<-(-log(-temp$mids[temp$mids<0],2)+log(0.00025,2))
  # this creates two peaks at zero, need to offset centers
  #xvals[1000]<-(-log)
  #xvals[1001]<-0.1
  plot(x = xvals, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)") 
  points(x = xvals, y = 1+log(temp.Fibroblasts.curated$counts,2), type = "h", col = "blue", lwd = 2.5)
  points(x = xvals, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "#FF000050", lwd = 2.5)
  }
dev.off()

# try "new" scaling method, using asinh(x) = ln(x + sqrt(x^2 + 1))
pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_asinh.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(asinh(POE.matrix[i,]), breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.pluripotents<-hist(asinh(POE.matrix[i,pluripotents]), breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.Fibroblasts.curated<-hist(asinh(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = seq(-1, 1, 0.001), plot = FALSE)

  plot(x = temp$mids, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)") 
  points(x = temp$mids, y = 1+log(temp.Fibroblasts.curated$counts,2), type = "h", col = "blue", lwd = 2.5)
  points(x = temp$mids, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "#FF000050", lwd = 2.5)
  }
dev.off()

which(POE.matrix[11,Fibroblasts.curated] == max(POE.matrix[11,Fibroblasts.curated], na.rm = TRUE))

# try signedlog scale with appropriate breaks
signedlog2<-function(x){sign(x)*log((1+1000*abs(x)), 2)}

pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_signedlog.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(signedlog2(POE.matrix[i,]), breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.pluripotents<-hist(signedlog2(POE.matrix[i,pluripotents]), breaks = seq(-1, 1, 0.001), plot = FALSE)
  temp.Fibroblasts.curated<-hist(signedlog2(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = seq(-1, 1, 0.001), plot = FALSE)

  plot(x = temp$mids, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)") 
  points(x = temp$mids, y = 1+log(temp.Fibroblasts.curated$counts,2), type = "h", col = "blue", lwd = 2.5)
  points(x = temp$mids, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "#FF000050", lwd = 2.5)
  }
dev.off()

signsq<-function(x){sign(x)*(x^2)}

glioblastoma<-unique(c(gsm.source[grep("glioblastoma", gsm.source[,2], ignore.case = TRUE),1],
                    gsm.source[grep("glioblastomas", gsm.source[,2], ignore.case = TRUE),1],
                    gsm.char[grep("glioblastoma", gsm.char[,2], ignore.case = TRUE),1],
                    gsm.char[grep("glioblastomas", gsm.char[,2], ignore.case = TRUE),1]
                              ))


neurosphere<-unique(c(gsm.source[grep("neurosphere", gsm.source[,2], ignore.case = TRUE),1],
                    gsm.source[grep("neurosphere", gsm.source[,2], ignore.case = TRUE),1],
                    gsm.char[grep("neurosphere", gsm.char[,2], ignore.case = TRUE),1],
                    gsm.char[grep("neurosphere", gsm.char[,2], ignore.case = TRUE),1]
                              ))


# possible expansion of pluripotent set
stemCellTypes<-read.delim("/data/shared/Fingerprint/GEO_ES_iPS_non_stem.txt", stringsAsFactors = FALSE)
pluripotents.2<-stemCellTypes[(stemCellTypes$SimpleCellType == "ES" |stemCellTypes$SimpleCellType == "iPS"),]
# remove MGU430B and HGU133B as they give strange profiles
pluripotents.2<-pluripotents.2$GSM[!(pluripotents.2$Platform == "GPL97" |pluripotents.2$Platform == "GPL340")]
pluripotents.2<-pluripotents.2[pluripotents.2 %in% colnames(POE.matrix)]
pluripotents.full<-union(pluripotents, pluripotents.2)

pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_nolog.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(signedlog2(POE.matrix[i,]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.pluripotents<-hist(signedlog2(POE.matrix[i,pluripotents.full]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.Fibroblasts.curated<-hist(signedlog2(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.glioblastoma<-hist(signedlog2(POE.matrix[i,glioblastoma[glioblastoma %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.neurosphere<-hist(signedlog2(POE.matrix[i,neurosphere[neurosphere %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  plot(x = temp$mids, y = 1+log(temp$counts,2), type = "h", main = rownames(POE.matrix)[i], xlab = "Geneset POE", col = "grey", lwd = 2.5, ylab = "1+log(frequency)", ylim = c(0,16), xaxt = "n")
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  points(x = temp$mids, y = 1+log(temp.Fibroblasts.curated$counts,2), type = "h", col = "#0000FF80", lwd = 2.5)
  points(x = temp$mids, y = 1+log(temp.pluripotents$counts,2), type = "h", col = "#FF000080", lwd = 2.5)
#  points(x = temp$mids, y = 1+log(temp.glioblastoma$counts,2), type = "h", col = "#00FF0050", lwd = 2.5)
#  points(x = temp$mids, y = 1+log(temp.neurosphere$counts,2), type = "h", col = "#00FFFF50", lwd = 2.5)
  }
dev.off()

pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_noylog.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(signedlog2(POE.matrix[i,]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.pluripotents<-hist(signedlog2(POE.matrix[i,pluripotents.full]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.Fibroblasts.curated<-hist(signedlog2(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.glioblastoma<-hist(signedlog2(POE.matrix[i,glioblastoma[glioblastoma %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.neurosphere<-hist(signedlog2(POE.matrix[i,neurosphere[neurosphere %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  par(mfcol = c(2,1))
  plot(x = temp$mids, y = temp$counts, type = "h", main = rownames(POE.matrix)[i], xlab = "", col = "grey", lwd = 2.5, ylab = "frequency", xaxt = "n")
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  plot(x = temp$mids, y = temp.Fibroblasts.curated$counts, type = "h", col = "#0000FF80", lwd = 2.5, xlab= "", xaxt = "n", ylab = "frequency", ylim = c(0, max(c(temp.Fibroblasts.curated$counts, temp.pluripotents$counts))))
  points(x = temp$mids, y = temp.pluripotents$counts, type = "h", col = "#FF000080", lwd = 2.5)
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
#  points(x = temp$mids, y = 1+log(temp.glioblastoma$counts,2), type = "h", col = "#00FF0050", lwd = 2.5)
#  points(x = temp$mids, y = 1+log(temp.neurosphere$counts,2), type = "h", col = "#00FFFF50", lwd = 2.5)
  }
dev.off()


# need to recurate the ES data to clean up some of the dubious iPS lines
# just create a new file ES/iPS.clean or something like that


# create lists of pluripotent arrays and fibroblast arrays used
wd<-getwd()
setwd("/home/galtschu2/Documents/Databases/GEOmetadb")
getSQLiteFile()
con <- dbConnect(SQLite(), "GEOmetadb.sqlite")
fibroblasts.curated.gse<-geoConvert(Fibroblasts.curated, c("gse"))$gse
fibroblasts.curated.gpl<-geoConvert(Fibroblasts.curated, c("gpl"))$gpl
Fibroblasts.curated.frame<-data.frame(GSM = Fibroblasts.curated,
                                      GSE = fibroblasts.curated.gse[match(Fibroblasts.curated, fibroblasts.curated.gse[,1]),2],
                                      GPL = fibroblasts.curated.gpl[match(Fibroblasts.curated, fibroblasts.curated.gpl[,1]),2],
                                      source = gsm.source[match(Fibroblasts.curated, gsm.source[,1]),2],
                                      Characteristics = gsm.char[match(Fibroblasts.curated, gsm.char[,1]),2]                                      
                                      )


pluripotents.full.gse<-geoConvert(pluripotents.full, c("gse"))$gse
pluripotents.full.gpl<-geoConvert(pluripotents.full, c("gpl"))$gpl
pluripotents.full.frame<-data.frame(GSM = pluripotents.full,
                                      GSE = pluripotents.full.gse[match(pluripotents.full, pluripotents.full.gse[,1]),2],
                                      GPL = pluripotents.full.gpl[match(pluripotents.full, pluripotents.full.gpl[,1]),2],
                                      source = gsm.source[match(pluripotents.full, gsm.source[,1]),2],
                                      Characteristics = gsm.char[match(pluripotents.full, gsm.char[,1]),2]
                                      )

dbDisconnect(con)
setwd(wd)
pluripotents.frame<-pluripotents.full.frame
fibroblasts.frame<-Fibroblasts.curated.frame
save(pluripotents.frame, file = "/data/shared/Fingerprint/curatedCellTypes/curatedPluripotentArrays.RData")
save(fibroblasts.frame, file = "/data/shared/Fingerprint/curatedCellTypes/curatedFibroblastsArrays.RData")

write.table(pluripotents.frame, file = "/data/shared/Fingerprint/curatedCellTypes/curatedPluripotentArrays.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(fibroblasts.frame, file = "/data/shared/Fingerprint/curatedCellTypes/curatedFibroblastArrays.txt",
sep = "\t", row.names = FALSE, quote = FALSE)


frame.curated<-as.data.frame(matrix(ncol = 1 + nrow(POE.matrix), nrow = length(c(rep("Fibroblast", length(Fibroblasts.curated)),
  rep("ES and iPS", length(pluripotents.full))
  ))))
colnames(frame.curated)<- c("Type", rownames(POE.matrix))
frame.curated$Type<-c(rep("Fibroblasts", length(Fibroblasts.curated)),
                      rep("ES and iPS", length(pluripotents.full))
                      )

# plot beanplot of the curated pluripotent and fibroblast arrays
for (i in 1:nrow(POE.matrix)){
  frame.curated[,(i+1)]<-c(
  signedlog2(POE.matrix[i, Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]),
  signedlog2(POE.matrix[i, pluripotents.full[pluripotents.full %in% colnames(POE.matrix)]])
  )
  }
  
pdf("pluriVsFibBean.pdf")
for (i in 1:nrow(POE.matrix)){
  beanplot(frame.curated[,i+1]~frame.curated[,1], las = 2, main = rownames(POE.matrix)[i]
          ,yaxt = "n", ylim = signedlog2(c(-1,1)), las = 1, what = c(FALSE, TRUE, TRUE, TRUE), beanline = "median")
  axis(2, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  
  }
dev.off()

# what is the mean of the ES cells
1-sum(POE.matrix[11,] > mean(POE.matrix[11,pluripotents.full]), na.rm = TRUE)/length(POE.matrix[11,])

which(frame.curated[,9+1] == max(frame.curated[,9+1]))

# replot using combined scale - somehow got the wrong matrix first time round..?
load("/data/shared/Fingerprint/Experimental.POE.matrix.2011-04-15")
# load re-curated list of pluripotent arrays
load("/data/shared/Fingerprint/curatedCellTypes/curatedPluripotentArrays_20110510.RData")

pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_plotstyle2.pdf")
for (i in 1:nrow(POE.matrix)){
  temp<-hist(signedlog2(POE.matrix[i,]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.pluripotents<-hist(signedlog2(POE.matrix[i,pluripotents.frame$GSM]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.Fibroblasts.curated<-hist(signedlog2(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.glioblastoma<-hist(signedlog2(POE.matrix[i,glioblastoma[glioblastoma %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.neurosphere<-hist(signedlog2(POE.matrix[i,neurosphere[neurosphere %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  # plot GEO corpus
  par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
  plot(x = temp$mids, y = temp$counts, type = "h", main = rownames(POE.matrix)[i], xlab = "", col = "grey", lwd = 2.5, ylab = "frequency", xaxt = "n", mar = c(0, 4, 4, 0))
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  # plot ES fibroblasts
  par(mar = c(7, 4, 4, 2))
  plot(x = temp$mids, y = temp.Fibroblasts.curated$counts, type = "h", col = "#0000FF80", lwd = 2.5, xlab= "Gene set expression score", xaxt = "n", ylab = "frequency", ylim = c(0, max(c(temp.Fibroblasts.curated$counts, temp.pluripotents$counts))))
  points(x = temp$mids, y = temp.pluripotents$counts, type = "h", col = "#FF000080", lwd = 2.5)
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
#  points(x = temp$mids, y = 1+log(temp.glioblastoma$counts,2), type = "h", col = "#00FF0050", lwd = 2.5)
#  points(x = temp$mids, y = 1+log(temp.neurosphere$counts,2), type = "h", col = "#00FFFF50", lwd = 2.5)
  }
dev.off()


# now just produce glioma.G, pathway 11
pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_plotstyle2_glioma.G.pdf")
for (i in 11){
  temp<-hist(signedlog2(POE.matrix[i,]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.pluripotents<-hist(signedlog2(POE.matrix[i,pluripotents.frame$GSM]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.Fibroblasts.curated<-hist(signedlog2(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.glioblastoma<-hist(signedlog2(POE.matrix[i,glioblastoma[glioblastoma %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.neurosphere<-hist(signedlog2(POE.matrix[i,neurosphere[neurosphere %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  # plot GEO corpus
  par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
  plot(x = temp$mids, y = temp$counts, type = "h", main = rownames(POE.matrix)[i], xlab = "", col = "grey", lwd = 2.5, ylab = "frequency", xaxt = "n", mar = c(0, 4, 4, 0))
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  # plot ES fibroblasts
  par(mar = c(7, 4, 4, 2))
  plot(x = temp$mids, y = temp.Fibroblasts.curated$counts, type = "h", col = "#0000FF80", lwd = 2.5, xlab= "Gene set expression score", xaxt = "n", ylab = "frequency", ylim = c(0, max(c(temp.Fibroblasts.curated$counts, temp.pluripotents$counts))))
  points(x = temp$mids, y = temp.pluripotents$counts, type = "h", col = "#FF000080", lwd = 2.5)
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
#  points(x = temp$mids, y = 1+log(temp.glioblastoma$counts,2), type = "h", col = "#00FF0050", lwd = 2.5)
#  points(x = temp$mids, y = 1+log(temp.neurosphere$counts,2), type = "h", col = "#00FFFF50", lwd = 2.5)
  }
dev.off()

# now just produce glioma.H, pathway 12
pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_plotstyle2_glioma.H.pdf")
for (i in 12){
  temp<-hist(signedlog2(POE.matrix[i,]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.pluripotents<-hist(signedlog2(POE.matrix[i,pluripotents.frame$GSM]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.Fibroblasts.curated<-hist(signedlog2(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.glioblastoma<-hist(signedlog2(POE.matrix[i,glioblastoma[glioblastoma %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
#  temp.neurosphere<-hist(signedlog2(POE.matrix[i,neurosphere[neurosphere %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  # plot GEO corpus
  par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
  plot(x = temp$mids, y = temp$counts, type = "h", main = rownames(POE.matrix)[i], xlab = "", col = "grey", lwd = 2.5, ylab = "frequency", xaxt = "n", mar = c(0, 4, 4, 0))
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  # plot ES fibroblasts
  par(mar = c(7, 4, 4, 2))
  plot(x = temp$mids, y = temp.Fibroblasts.curated$counts, type = "h", col = "#0000FF80", lwd = 2.5, xlab= "Gene set expression score", xaxt = "n", ylab = "frequency", ylim = c(0, max(c(temp.Fibroblasts.curated$counts, temp.pluripotents$counts))))
  points(x = temp$mids, y = temp.pluripotents$counts, type = "h", col = "#FF000080", lwd = 2.5)
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
#  points(x = temp$mids, y = 1+log(temp.glioblastoma$counts,2), type = "h", col = "#00FF0050", lwd = 2.5)
#  points(x = temp$mids, y = 1+log(temp.neurosphere$counts,2), type = "h", col = "#00FFFF50", lwd = 2.5)
  }
dev.off()

# calculate empirically derived p-values

sum(POE.matrix[11,] > median(POE.matrix[11,pluripotents.frame$GSM]), na.rm = TRUE)/length(na.omit(POE.matrix[11,]))
sum(POE.matrix[11,] > median(POE.matrix[12,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), na.rm = TRUE)/length(na.omit(POE.matrix[12,]))
sum(POE.matrix[12,] > median(POE.matrix[12,pluripotents.frame$GSM]), na.rm = TRUE)/length(na.omit(POE.matrix[12,]))
sum(POE.matrix[12,] > median(POE.matrix[12,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), na.rm = TRUE)/length(na.omit(POE.matrix[12,]))

#######
# 25th June 2011 - Cecilie edits
# Repeat analysis, separating out the ES and iPS signature
pluripotents.frame$GSM<-as.character(pluripotents.frame$GSM)
ES.terms<-c("embryonic", "ES", "ESC", "Fetal")
iPS.terms<-c("iPS", "iPSC", "induced", "MCV8.1")

iPS.frame<-pluripotents.frame[(rowSums(sapply(iPS.terms, function(x){
    grepl(x, apply(pluripotents.frame[,c("Characteristics", "source")], 1, paste, collapse = " "), ignore.case = TRUE)
    })
    ) > 0),]
ES.frame<-pluripotents.frame[(rowSums(sapply(ES.terms, function(x){
    grepl(x, apply(pluripotents.frame[,c("Characteristics", "source")], 1, paste, collapse = " "), ignore.case = TRUE)
    })
    ) > 0),]

# which are missed?
pluripotents.frame[match(setdiff(as.character(pluripotents.frame$GSM), union(iPS.frame$GSM, ES.frame$GSM)), pluripotents.frame$GSM),1:5]
# which are in both?
pluripotents.frame[match(intersect(iPS.frame$GSM, ES.frame$GSM), pluripotents.frame$GSM),1:5]
# these should all only be in the iPS panel
ES.frame<-ES.frame[-match(intersect(iPS.frame$GSM, ES.frame$GSM), ES.frame$GSM),]

# okay, now re-plot histograms
pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_plotstyle2_glioma.G_iPS_ES.pdf")
for (i in 11){
  temp<-hist(signedlog2(POE.matrix[i,]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.ES<-hist(signedlog2(POE.matrix[i,ES.frame$GSM]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.iPS<-hist(signedlog2(POE.matrix[i,iPS.frame$GSM]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.Fibroblasts.curated<-hist(signedlog2(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  # plot GEO corpus
  par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
  plot(x = temp$mids, y = temp$counts, type = "h", main = rownames(POE.matrix)[i], xlab = "", col = "grey", lwd = 2.5, ylab = "frequency", xaxt = "n", mar = c(0, 4, 4, 0))
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  # plot ES fibroblasts
  par(mar = c(7, 4, 4, 2))
  plot(x = temp$mids, y = temp.Fibroblasts.curated$counts, type = "h", col = "#0000FF80", lwd = 2.5, xlab= "Gene set expression score", xaxt = "n", ylab = "frequency", ylim = c(0, max(c(temp.Fibroblasts.curated$counts, temp.pluripotents$counts))))
  points(x = temp$mids, y = temp.iPS$counts, type = "h", col = "#00FF0080", lwd = 2.5)
  points(x = temp$mids, y = temp.ES$counts, type = "h", col = "#FF000080", lwd = 2.5)
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  }
dev.off()

pdf("Experimental_POE_histograms_with_pluripotents_and_Fibroblasts.curated_plotstyle2_glioma.H_iPS_ES.pdf")
for (i in 12){
  temp<-hist(signedlog2(POE.matrix[i,]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.ES<-hist(signedlog2(POE.matrix[i,ES.frame$GSM]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.iPS<-hist(signedlog2(POE.matrix[i,iPS.frame$GSM]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  temp.Fibroblasts.curated<-hist(signedlog2(POE.matrix[i,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), breaks = (seq(-10, 10, 0.1)), plot = FALSE)
  # plot GEO corpus
  par(mfcol = c(2,1), mar = c(0, 4, 4, 2))
  plot(x = temp$mids, y = temp$counts, type = "h", main = rownames(POE.matrix)[i], xlab = "", col = "grey", lwd = 2.5, ylab = "frequency", xaxt = "n", mar = c(0, 4, 4, 0))
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  # plot ES fibroblasts
  par(mar = c(7, 4, 4, 2))
  plot(x = temp$mids, y = temp.Fibroblasts.curated$counts, type = "h", col = "#0000FF80", lwd = 2.5, xlab= "Gene set expression score", xaxt = "n", ylab = "frequency", ylim = c(0, max(c(temp.Fibroblasts.curated$counts, temp.pluripotents$counts))))
  points(x = temp$mids, y = temp.iPS$counts, type = "h", col = "#00FF0080", lwd = 2.5)
  points(x = temp$mids, y = temp.ES$counts, type = "h", col = "#FF000080", lwd = 2.5)
  axis(1, at = signedlog2(c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1)), labels = c(-1,-0.1, -0.01, 0, 0.01, 0.1, 1))
  }
dev.off()


sum(POE.matrix[11,] > median(POE.matrix[11,ES.frame$GSM]), na.rm = TRUE)/length(na.omit(POE.matrix[11,]))
sum(POE.matrix[11,] > median(POE.matrix[11,iPS.frame$GSM]), na.rm = TRUE)/length(na.omit(POE.matrix[11,]))
sum(POE.matrix[11,] > median(POE.matrix[12,Fibroblasts.curated[Fibroblasts.curated %in% colnames(POE.matrix)]]), na.rm = TRUE)/length(na.omit(POE.matrix[12,]))

write.table(ES.frame, file = "/data/shared/Fingerprint/curatedCellTypes/curatedESArrays.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(iPS.frame, file = "/data/shared/Fingerprint/curatedCellTypes/curatediPSArrays.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

save.image("Experimental.POE.Analysis.RData")
