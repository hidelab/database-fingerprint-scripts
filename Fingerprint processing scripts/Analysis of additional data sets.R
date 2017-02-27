# Analysis of new experimental data in fingerprint system

# Strategy
# 1) Normalize fingerprint as before
# 2) Sample arrays at stringent threshold level (5%)
# 3) Correlate with other pathways

# Local

setwd("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Data")
load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Data/analysisOfAdditionalDataSets.RData")

source("/Users/GabrielAltschuler/Dropbox/gabriels-scripts/gabriel functions.R")


# Server loading and analysis of fingerprints
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/")

fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints"
setwd(fingerpath)
source("/home/galtschu2/Documents/Databases/gabriel functions.R")

if (!(exists("ExperimentalGeneSigs.Hs.gs"))){
	load("/home/galtschu2/Documents/Databases/Gene sets/ExperimentalGeneSigs.Hs.gs", .GlobalEnv)
	}
if (!(exists("ExperimentalGeneSigs.Mm.gs"))){
	load("/home/galtschu2/Documents/Databases/Gene sets/ExperimentalGeneSigs.Mm.gs", .GlobalEnv)
	}

if (!(exists("chipframe"))){
	load("/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update.RData", .GlobalEnv)
	}


# only load complete fingerprints this takes some time (10mins) on the server
# Incomplete run - currently have 93994 files

files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
files<-files[grep("Experimental", files)]

frame<-vector("list", length(files))
pb <- txtProgressBar(min = 0, max = length((frame)), style = 3)
for (i in 1:length(files)){
	load(files[i])
	setTxtProgressBar(pb, i)
	frame[i]<-temp1
	}

SCG<-vector("list", length(frame))
pb <- txtProgressBar(min = 0, max = length((frame)), style = 3)
for (i in 1:length(frame)){
	SCG[i]<-frame[[i]]$SCG
	setTxtProgressBar(pb, i)
	}
	
names(SCG)<-gsub("Experimental_", "", files)
SCG.frame<-t(as.data.frame(SCG))
colnames(SCG.frame)<-names(ExperimentalGeneSigs.Hs.gs)

platform<-vector("list", length(frame))
pb <- txtProgressBar(min = 0, max = length((frame)), style = 3)
for (i in 1:length(frame)){
	platform[i]<-frame[[i]]$platform
	setTxtProgressBar(pb, i)
	}
	
chiplengths <- c(16477, 12170, 5779, 18110, 12056, 7738, 8341, 12170, 12056, 7030, length(unique(chipframe$GPL2986$ann$EntrezID)))
names(chiplengths)<-c("GPL1261", "GPL339", "GPL340", "GPL570", "GPL571", "GPL81", "GPL8300", "GPL8321", "GPL96", "GPL97", "GPL2986")
species <- c("mouse", "mouse", "mouse", "human", "human", "mouse", "human", "mouse", "human", "human", "human")
names(species)<-names(chiplengths)

platform.lengths <- vector("list", length(platform))
pb <- txtProgressBar(min = 0, max = length((frame)), style = 3)
for (i in 1:length(platform.lengths)){
	platform.lengths[i]<-chiplengths[platform[[i]]]
	setTxtProgressBar(pb, i)
	}

chipgenes<-vector("list", length(chiplengths))
names(chipgenes)<-names(chiplengths)
for (i in 1:length(chipgenes)){
	chipgenes[i]<-list(unique(gsub("_at", "", chipframe[[names(chiplengths)[i]]]$ann$EntrezID)))
	}


ExperimentalGeneSigs.Hs.overlaplength<-function(x){
	temp<-vector("numeric", length(ExperimentalGeneSigs.Hs.gs))
	for (i in 1:length(ExperimentalGeneSigs.Hs.gs)){
		temp[i]<-length(intersect(x, ExperimentalGeneSigs.Hs.gs[[i]]))		}
	return(temp)
	}

ExperimentalGeneSigs.Mm.overlaplength<-function(x){
	temp<-vector("numeric", length(ExperimentalGeneSigs.Mm.gs))
	for (i in 1:length(ExperimentalGeneSigs.Mm.gs)){
		temp[i]<-length(intersect(x, ExperimentalGeneSigs.Mm.gs[[i]]))		}
	return(temp)
	}

overlaps<-vector("list", length(chiplengths))
names(overlaps)<-names(chiplengths)
for (i in 1:length(chipgenes)){
	if (species[names(overlaps[i])] == "human"){
		overlaps[i]<-list(ExperimentalGeneSigs.Hs.overlaplength(chipgenes[[names(chiplengths)[i]]]))
		}
	if (species[names(overlaps[i])] == "mouse"){
		overlaps[i]<-list(ExperimentalGeneSigs.Mm.overlaplength(chipgenes[[names(chiplengths)[i]]]))
		} 
	}

# first use the model to fit the mean and the sd for each chip

mean.m <- 0.0555961
mean.c <- 0.1135812
sd.b.m <- 0.2734675
sd.b.c <- 58.91614
par.sd.a <- 9.812654
sd.power.m <- (1.042497e-06)
sd.power.c <- (-0.5250856)

par.mean <- as.numeric((mean.m*chiplengths)+mean.c)
par.sd.b<-as.numeric((sd.b.m*chiplengths)+sd.b.c)
par.sd.power<-as.numeric((sd.power.m*chiplengths)+sd.power.c)

names(par.mean) <- names(chiplengths)
names(par.sd.b) <- names(chiplengths)
names(par.sd.power) <- names(chiplengths)

# Don't want NA's for now so relax threshold

Norm.SCG<-vector("list", length(SCG))
pb <- txtProgressBar(min = 0, max = length(SCG), style = 3)
for (i in 1:length(SCG)){
	intersect <- overlaps[[platform[[i]]]]
	intersect[intersect<3] <- NA
	mean<-par.mean[platform[[i]]]*(1-(intersect^-1))
	sd<- par.sd.a+(par.sd.b[platform[[i]]]*(intersect^par.sd.power[platform[[i]]]))
	Norm.SCG[i]<-list(pnorm(SCG[[i]], mean = mean, sd = sd))
	setTxtProgressBar(pb, i)
	}

names(Norm.SCG)<-gsub("Experimental_", "", files)

# problem with zero length arrays
Norm.SCG.lengths<-lapply(Norm.SCG, length)
zeros<-names(Norm.SCG[unlist(Norm.SCG.lengths) == 0])
Norm.SCG[match(zeros, names(Norm.SCG))]<-list(rep(1, length(ExperimentalGeneSigs.Hs.gs)))

Norm.SCG.frame<-t(as.data.frame(Norm.SCG))
colnames(Norm.SCG.frame)<-names(ExperimentalGeneSigs.Hs.gs)
save(Norm.SCG.frame, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/Norm.SCG.experimental.frame.RData")

# normalize per chip
Norm.SCG.matrix<-as.matrix(t(Norm.SCG.frame))
for (i in 1:length(chiplengths)){
	temp<-Norm.SCG.matrix[,grep(names(chiplengths)[i], platform)]
	assign(paste(names(chiplengths)[i], ".Norm.SCG.matrix", sep = ""), temp)
	}
	
for (i in 1:length(chiplengths)){
	temp<-apply(get(paste(names(chiplengths)[i], ".Norm.SCG.matrix", sep = "")), 1, ecdf)
	assign(paste(names(chiplengths)[i], ".Norm.SCG.ecdf", sep = ""), temp)
	save(temp, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", names(chiplengths)[i], ".experimental_Norm.SCG.ecdf.RData", sep = ""))
	}


breaks<-0.005*(0:200)
pdf(file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/Norm_hist_by_chip_experimental.pdf")
par(mfcol = c(4,2))
for (j in 1:length(Norm.SCG.matrix[,1])){
	for (i in c(1,2,4,5,6,7,8,9)){	
		hist(get(paste(names(chiplengths)[i], ".Norm.SCG.matrix", sep = ""))[j,], xlim = c(0,1), breaks = breaks, main = paste(names(chiplengths)[i], rownames(Norm.SCG.matrix)[j], sep = "_"), xlab = "Q1 normalized score")
		}
	}
	
dev.off()

# now normalize each chip separately
for (i in 1:length(chiplengths)){
	temp1<-get(paste(names(chiplengths)[i], ".Norm.SCG.matrix", sep = ""))
	temp<- matrix(ncol = length(temp1[1,]), nrow = length(temp1[,1]))
	colnames(temp)<-colnames(temp1)
	rownames(temp)<-rownames(temp1)
	temp2<-get(paste(names(chiplengths)[i], ".Norm.SCG.ecdf", sep = ""))
	for (j in 1:length(temp2)){
		temp[j,]<-temp2[[j]](temp1[j,])
		}				
	assign(paste(names(chiplengths)[i], ".N.Norm.SCG", sep = ""), temp)
	}

# re-plotting should all be uniform distributions (except the NAs) - it worked
pdf(file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/Norm_Norm_hist_by_chip_experimental.pdf")	
par(mfcol = c(4,2))
for (j in 1:length(Norm.SCG.matrix[,1])){
	for (i in c(1,2,4,5,6,7,8,9)){	
		hist(get(paste(names(chiplengths)[i], ".N.Norm.SCG", sep = ""))[j,], xlim = c(0,1), breaks = breaks, main = paste(names(chiplengths)[i], rownames(Norm.SCG.matrix)[j], sep = "_"), xlab = "Q1 normalized score")
		}
	}
	
dev.off()

# re-assemble matrix_do not include ABI for now

ReNorm.matrix<-cbind(
					get(paste(names(chiplengths)[1], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[2], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[3], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[4], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[5], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[6], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[7], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[8], ".N.Norm.SCG", sep = "")),					get(paste(names(chiplengths)[9], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[10], ".N.Norm.SCG", sep = ""))
					
					)

#save(ReNorm.matrix, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix_experimental.RData")
save(ReNorm.matrix, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix_experimental_full.RData")
##################
Now on local machine
load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/GEOdescription.RData")
setwd("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Analysis")
load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/data/Norm.SCG.experimental.frame.RData")
load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/data/ReNorm.matrix_experimental.RData")
load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/data/platform.frame.RData")
description[match(colnames(ReNorm.matrix[, ReNorm.matrix[1,] < 0.02 & ReNorm.matrix[2,] > 0.98]), names(description))]
platform.frame[match(colnames(ReNorm.matrix[, ReNorm.matrix[1,] < 0.05 & ReNorm.matrix[2,] > 0.95]), rownames(platform.frame)),]

description[match(colnames(ReNorm.matrix[,ReNorm.matrix[2,] > 0.999]), names(description))]
platform.frame[match(colnames(ReNorm.matrix[,ReNorm.matrix[2,] > 0.999]), rownames(platform.frame)),]

# where do pluripotents sit..?

setwd("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Data")
database<-read.delim("GEO_ES_iPS_non_stem.txt")
database.valid<-database[(database$GSM %in% colnames(ReNorm.matrix)),]
database.human<-database.valid[grep("Human", database.valid$Species),]
# remove outlier arrays
outliers<-c("GSM423940", "GSM423941", "GSM423942", "GSM423943", "GSM423944", "GSM423945")
database.human<-database.human[!(database.human$GSM %in% outliers),]
database.valid<-database.valid[!(database.valid$GSM %in% outliers),]
pluripotents<-as.character(database.human[c(grep("ES", database.human$SimpleCellType) ,grep("^iPS", database.human$SimpleCellType)),3])
pluripotents.all<-as.character(database.valid[c(grep("ES", database.valid$SimpleCellType) ,grep("^iPS", database.valid$SimpleCellType)),3])
par(mfcol = c(2,1))
hist(ReNorm.matrix[1,pluripotents], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[1], xlab = "normalized score")
hist(ReNorm.matrix[2,pluripotents], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[2], xlab = "normalized score")
par(mfcol = c(2,1))
hist(ReNorm.matrix[1,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[1], xlab = "normalized score")
hist(ReNorm.matrix[2,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[2], xlab = "normalized score")

# what about the Ovarian arrays themselves..?

ovarian.list<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Ovarian cancer study/Ovarian cancer GSM list.txt")
ovarians<-as.character(ovarian.list$GEO)
hist(ReNorm.matrix[1,ovarians[ovarians %in% colnames(ReNorm.matrix)]], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[1], xlab = "normalized score")
hist(ReNorm.matrix[2,ovarians[ovarians %in% colnames(ReNorm.matrix)]], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[2], xlab = "normalized score")
plot(ReNorm.matrix[1,ovarians[ovarians %in% colnames(ReNorm.matrix)]], ReNorm.matrix[2,ovarians[ovarians %in% colnames(ReNorm.matrix)]], xlab = rownames(ReNorm.matrix)[1], ylab = rownames(ReNorm.matrix)[2])

# histograms for other gene sets - CD133 normalization did not work
hist(ReNorm.matrix[3,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[3], xlab = "normalized score")
hist(ReNorm.matrix[4,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[4], xlab = "normalized score")

# Interestingly the breast cancer siRNA hits are upregulated in ES/iPS cells
hist(ReNorm.matrix[5,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[5], xlab = "normalized score")
hist(ReNorm.matrix[6,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[6], xlab = "normalized score")

# Pluripotent cells show low expression in Cotton Workers gene list
hist(ReNorm.matrix[7,], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[7], xlab = "normalized score")
hist(ReNorm.matrix[7,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[7], xlab = "normalized score")

# Not much change in enrichment for miR34a list
hist(ReNorm.matrix[8,], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[8], xlab = "normalized score")
hist(ReNorm.matrix[8,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[8], xlab = "normalized score")

# Striking pattern for GBM gene lists
hist(ReNorm.matrix[9,], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[9], xlab = "normalized score")
hist(ReNorm.matrix[10,], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[10], xlab = "normalized score")

hist(ReNorm.matrix[9,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[9], xlab = "normalized score")
hist(ReNorm.matrix[10,pluripotents.all], xlim = c(0,1), breaks  =seq(0,1,0.01), main = rownames(ReNorm.matrix)[10], xlab = "normalized score")

# what are the arrays that show highest expression of the Langmoen GBM upregualted signature? 875 arrays have score = 1
GBM.high<-colnames(ReNorm.matrix)[ReNorm.matrix[9,] == 1]
GBM.high_GBM.low<-colnames(ReNorm.matrix)[ReNorm.matrix[9,] == 1 & ReNorm.matrix[10,] < 0.1]

load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/GEOdescription.RData")

GBM.frame<-data.frame(GEO = GBM.high_GBM.low, description = unlist(description[GBM.high_GBM.low]))

# what is the correlation in the fingerprints of the signatures..?
# do not include the CD133 as they didn't normalize properly

heatmap(cor(t(ReNorm.matrix)), margins = c(10,10), cexRow = 0.75, cexCol = 0.75)

heatmap.2(cor(t(ReNorm.matrix[-c(3,4),])), margins = c(15,15), cexRow = 0.75, cexCol = 0.75, col = colorRampPalette(colors = c("blue", "grey", "red"), bias = 1.35), trace = "none", density.info = "none")

heatmap.2(cor(t(ReNorm.matrix[, pluripotents.all])), margins = c(10,10), cexRow = 0.75, cexCol = 0.75, col = colorRampPalette(c("blue", "grey", "red")), trace = "none", density.info = "none")

# how dies this fit-in with the other pathways..?

load("ReNorm.matrix.15ternary_7_10_10.RData")

col.intersect<-intersect(colnames(ReNorm.matrix.ternary.15), colnames(ReNorm.matrix))

fullmatrix<-rbind(ReNorm.matrix.ternary.15[, col.intersect], ReNorm.matrix[, col.intersect])

cor<-cor(t(fullmatrix))
zeropathways <-c(
	"00300 Lysine biosynthesis",
	"00401 Novobiocin biosynthesis",
	"00471 D-Glutamine and D-glutamate metabolism",
	"00472 D-Arginine and D-ornithine metabolism", 
	"00643 Styrene degradation",                  
	"00780 Biotin metabolism",                     
	"00785 Lipoic acid metabolism",                
	"amino acid conjugation of benzoic acid",      
	"Polyol pathway",                              
	"Serotonin Receptor 2 -> STAT3 signaling",     
	"Alpha6 Beta4 Integrin DOWN",                  
	"Hedgehog DOWN",                               
	"Kit Receptor DOWN",                           
	"IL-9 DOWN"
	)   

cor["LangmoenGBM_Down",abs(cor["LangmoenGBM_Down",])>0.4 & is.na(abs(cor["LangmoenGBM_Down",])) == FALSE]


GMB.down.cor<-rownames(cor)[abs(cor["LangmoenGBM_Down",])>0.45 & is.na(abs(cor["LangmoenGBM_Down",])) == FALSE]
GMB.up.cor<-rownames(cor)[abs(cor["LangmoenGBM_Up",])>0.45 & is.na(abs(cor["LangmoenGBM_Up",])) == FALSE]

GBM.cor<-unique(c(GMB.down.cor, GMB.up.cor))

library(gplots)
heatmap.2((cor[c("LangmoenGBM_Up", "LangmoenGBM_Down"), GBM.cor]), margins = c(20,15), scale = "none", col = colorRampPalette(colors = c("blue", "grey", "red"), bias = 1.35), Rowv = FALSE, cexCol = 0.7, cexRow = 0.7, trace = "none", density.info = "none")

# compare this to the number of overlapping genes

load("/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Hs.gs")

load("/Users/GabrielAltschuler/Documents/Databases/Gene sets/kegg_wiki_TR_static.Hs.gs")

kegg_wiki_TR_static.Hs.gs<-append(kegg_wiki_TR_static.Hs.gs, ExperimentalGeneSigs.Hs.gs[-c(9,10)])

GMB.down.intersects<-vector("numeric", length(kegg_wiki_TR_static.Hs.gs))
for (i in 1:length(kegg_wiki_TR_static.Hs.gs)){
	GMB.down.intersects[i]<-length(intersect(ExperimentalGeneSigs.Hs.gs[[10]], kegg_wiki_TR_static.Hs.gs[[i]]))
	}

GMB.down.intersects<-hyperPathway(ExperimentalGeneSigs.Hs.gs[[10]], kegg_wiki_TR_static.Hs.gs, 20000)

GMB.up.intersects<-hyperPathway(ExperimentalGeneSigs.Hs.gs[[9]], kegg_wiki_TR_static.Hs.gs, 20000)

GBM.intersects<-data.frame("GMB_up" = GMB.up.intersects$nGenes,"GMB_down" = GMB.down.intersects$nGenes)
rownames(GBM.intersects)<-names(kegg_wiki_TR_static.Hs.gs)			
GBM.pvals<-data.frame("GMB_up" = -log(GMB.up.intersects$"P-value",10),"GMB_down" = -log(GMB.down.intersects$"P-value", 10))
rownames(GBM.pvals)<-names(kegg_wiki_TR_static.Hs.gs)



temp<-heatmap.2((cor[c("LangmoenGBM_Up", "LangmoenGBM_Down"), GBM.cor]))

heatmap.2(t(as.matrix(GBM.intersects[rownames(temp$carpet),])), margins = c(20,15), scale = "none", col = colorRampPalette(colors = c("light grey", "red")), dendrogram = "none", Colv = FALSE, Rowv = FALSE, cexCol = 0.7, cexRow = 0.7, trace = "none", density.info = "none")

heatmap.2(t(as.matrix(GBM.pvals[rownames(temp$carpet),])), margins = c(20,15), scale = "none", col = colorRampPalette(colors = c("light grey", "red", "red", "red")), dendrogram = "none", Colv = FALSE, Rowv = FALSE, cexCol = 0.7, cexRow = 0.7, trace = "none", density.info = "none")

##############
# Comparison with other proliferative cell types for cotton workers
####### Background for ES cells

load("/Users/GabrielAltschuler/Documents/Databases/GEOmetadb/gsm_characteristics.RData")

HeLa<-gsm.char[grep("HeLa", gsm.char[,2]),1]
Jurkat<-gsm.char[grep("Jurkat", gsm.char[,2]),1]
MCF.7<-gsm.char[grep("MCF-7", gsm.char[,2]),1]
HEK.293<-gsm.char[grep("HEK293", gsm.char[,2]),1]
iPS<-gsm.char[grep("iPS", gsm.char[,2]),1]
U2OS<-gsm.char[grep("U2OS", gsm.char[,2]),1]
Caco.2<-gsm.char[grep("Caco-2", gsm.char[,2]),1]
HT29<-gsm.char[grep("HT29", gsm.char[,2]),1]
T.47D<-gsm.char[grep("T47D", gsm.char[,2]),1]
Mm.3T3<-gsm.char[grep("3T3", gsm.char[,2]),1]



prolif.frame.CottonWorkersMouseModel<-data.frame(
					MARS = c(
							ReNorm.matrix["CottonWorkersMouseModel", HeLa[HeLa %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", Jurkat[Jurkat %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", MCF.7[MCF.7 %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", HEK.293[HEK.293 %in% colnames(ReNorm.matrix)]],
						#	ReNorm.matrix["CottonWorkersMouseModel", iPS[iPS %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", U2OS[U2OS %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", Caco.2[Caco.2 %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", HT29[HT29 %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", pluripotents[pluripotents %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", Mm.3T3[Mm.3T3 %in% colnames(ReNorm.matrix)]],
							ReNorm.matrix["CottonWorkersMouseModel", T.47D[T.47D %in% colnames(ReNorm.matrix)]]
							),
					Type = c(
							rep(paste("HeLa (", sum(HeLa %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(HeLa %in% colnames(ReNorm.matrix))),
							rep(paste("Jurkat (", sum(Jurkat %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(Jurkat %in% colnames(ReNorm.matrix))),
							rep(paste("MCF.7 (", sum(MCF.7 %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(MCF.7 %in% colnames(ReNorm.matrix))),
							rep(paste("HEK.293 (", sum(HEK.293 %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(HEK.293 %in% colnames(ReNorm.matrix))),
						#	rep(paste("iPS (", sum(iPS %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(iPS %in% colnames(ReNorm.matrix))),
							rep(paste("U2OS (", sum(U2OS %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(U2OS %in% colnames(ReNorm.matrix))),
							rep(paste("Caco.2 (", sum(Caco.2 %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(Caco.2 %in% colnames(ReNorm.matrix))),
							rep(paste("HT29 (", sum(HT29 %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(HT29 %in% colnames(ReNorm.matrix))),
							rep(paste("ES/iPS (", sum(pluripotents %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(pluripotents %in% colnames(ReNorm.matrix))),
							rep(paste("3T3 (", sum(Mm.3T3 %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(Mm.3T3 %in% colnames(ReNorm.matrix))),
							rep(paste("T.47D (", sum(T.47D %in% colnames(ReNorm.matrix)), ")", sep = ""), sum(T.47D %in% colnames(ReNorm.matrix)))
							)
					)
pdf("/Users/GabrielAltschuler/Documents/Projects/CottonWorkers/CottonWorkersBoxplot.pdf")
par(mar = c(10,4,4,2))
boxplot(prolif.frame.CottonWorkersMouseModel[,1]~prolif.frame.CottonWorkersMouseModel[,2], las = 2, ylab = "Normalized CottonWorkersMouseModel expression", ylim = c(0,1))
dev.off()

# correlation of cotton workers with other pathways in the fingerprint
cor["CottonWorkersMouseModel",abs(cor["CottonWorkersMouseModel",])>0.45 & is.na(abs(cor["CottonWorkersMouseModel",])) == FALSE]
CottonWorkersMouseModel.cor<-rownames(cor)[abs(cor["CottonWorkersMouseModel",])>0.45 & is.na(abs(cor["CottonWorkersMouseModel",])) == FALSE]
# remove the cotten workers pathway itself
CottonWorkersMouseModel.cor<-CottonWorkersMouseModel.cor[-length(CottonWorkersMouseModel.cor)]
CottonWorkersMouseModel.hyper<-hyperPathway(ExperimentalGeneSigs.Hs.gs[["CottonWorkersMouseModel"]], kegg_wiki_TR_static.Hs.gs, 20000)

CottonWorkersMouseModel.intersects<-data.frame("CottonWorkersMouseModel" = CottonWorkersMouseModel.hyper$nGenes)
rownames(CottonWorkersMouseModel.intersects)<-names(kegg_wiki_TR_static.Hs.gs)			
CottonWorkersMouseModel.pvals<-data.frame("CottonWorkersMouseModel" = -log(CottonWorkersMouseModel.hyper$"P-value",10))
rownames(CottonWorkersMouseModel.pvals)<-names(kegg_wiki_TR_static.Hs.gs)

# plot heatmap
#heatmap.2(cbind(
#				cor["CottonWorkersMouseModel", CottonWorkersMouseModel.cor],
#				CottonWorkersMouseModel.intersects[CottonWorkersMouseModel.cor,],
#				CottonWorkersMouseModel.pvals[CottonWorkersMouseModel.cor,]
#				),
#				margins = c(20,15), scale = "col", col = colorRampPalette(colors = c("blue", "grey", "red"), bias = 1.35), Rowv = FALSE, cexCol = 0.7, cexRow = 0.7, trace = "none", density.info = "none")
# Doesn't really work to plot them all together, try separates but double up to allow heatmap
# Correlation
pdf("/Users/GabrielAltschuler/Documents/Projects/CottonWorkers/CottonWorkersCorrelationHeatmap.pdf")
CottonWorkersMouseModel.heatmap<-heatmap.2(cbind(
				cor["CottonWorkersMouseModel", CottonWorkersMouseModel.cor],
				cor["CottonWorkersMouseModel", CottonWorkersMouseModel.cor]
				),
				margins = c(5,29), scale = "none", col = colorRampPalette(colors = c("grey", "red")),
				cexCol = 0.7, cexRow = 0.7, trace = "none", density.info = "none"
				)
dev.off()
# Intersects
cotton.intersects<-cbind(
				CottonWorkersMouseModel.intersects[colnames(CottonWorkersMouseModel.heatmap$carpet),,drop = FALSE],
				CottonWorkersMouseModel.intersects[colnames(CottonWorkersMouseModel.heatmap$carpet),]
				)
pdf("/Users/GabrielAltschuler/Documents/Projects/CottonWorkers/CottonWorkersIntersectHeatmap.pdf")
heatmap.2(as.matrix(cotton.intersects),
				margins = c(5,29), scale = "none", col = colorRampPalette(colors = c("grey", "red")),
				cexCol = 0.7, cexRow = 0.7, trace = "none", density.info = "none",
				dendrogram = "none", Colv = FALSE, Rowv = FALSE
				)
dev.off()

# pvalues
cotton.pvals<-cbind(
				CottonWorkersMouseModel.pvals[colnames(CottonWorkersMouseModel.heatmap$carpet),,drop = FALSE],
				CottonWorkersMouseModel.pvals[colnames(CottonWorkersMouseModel.heatmap$carpet),]
				)
pdf("/Users/GabrielAltschuler/Documents/Projects/CottonWorkers/CottonWorkersPvalsHeatmap.pdf")
heatmap.2(as.matrix(cotton.pvals),
				margins = c(5,29), scale = "none", col = colorRampPalette(colors = c("grey", "red")),
				cexCol = 0.7, cexRow = 0.7, trace = "none", density.info = "none",
				dendrogram = "none", Colv = FALSE, Rowv = FALSE
				)
dev.off()


# what are the arrays that show highest expression of the cotton workers signature?
CottonWorkersMouseModel.high<-colnames(ReNorm.matrix)[ReNorm.matrix["CottonWorkersMouseModel",] == 1]
CottonWorkersMouseModel.high.desc<-gsm.char[match(CottonWorkersMouseModel.high, gsm.char$gsm),]
# how many of these contain the word lung?
length(union(grep("lung", CottonWorkersMouseModel.high.desc$characteristics_ch1), grep("Lung", CottonWorkersMouseModel.high.desc$characteristics_ch1)))


CottonWorkersMouseModel.high.string<-""
for (i in 1:length(unlist(CottonWorkersMouseModel.high.desc$characteristics_ch1))){
	CottonWorkersMouseModel.high.string<-paste(CottonWorkersMouseModel.high.string, unlist(CottonWorkersMouseModel.high.desc$characteristics_ch1)[i], sep = " ")
	}
CottonWorkersMouseModel.high.string<-gsub("[[:punct:]]", " ", CottonWorkersMouseModel.high.string)
CottonWorkersMouseModel.high.string<-gsub("\t", " ", CottonWorkersMouseModel.high.string)
CottonWorkersMouseModel.high.string<-toLower(CottonWorkersMouseModel.high.string)
CottonWorkersMouseModel.high.table<-table(strsplit(CottonWorkersMouseModel.high.string, " "))

head(sort(CottonWorkersMouseModel.high.table, decreasing = TRUE),60)

CottonWorkersMouseModel.high.95<-colnames(ReNorm.matrix)[ReNorm.matrix["CottonWorkersMouseModel",] > 0.95]
CottonWorkersMouseModel.high.95.desc<-gsm.char[match(CottonWorkersMouseModel.high.95, gsm.char$gsm),]
# how many of these contain the word lung?
length(union(grep("lung", CottonWorkersMouseModel.high.95.desc$characteristics_ch1), grep("Lung", CottonWorkersMouseModel.high.95.desc$characteristics_ch1)))


CottonWorkersMouseModel.high.95.string<-""
for (i in 1:length(unlist(CottonWorkersMouseModel.high.95.desc$characteristics_ch1))){
	CottonWorkersMouseModel.high.95.string<-paste(CottonWorkersMouseModel.high.95.string, unlist(CottonWorkersMouseModel.high.95.desc$characteristics_ch1)[i], sep = " ")
	}
CottonWorkersMouseModel.high.95.string<-gsub("[[:punct:]]", " ", CottonWorkersMouseModel.high.95.string)
CottonWorkersMouseModel.high.95.string<-gsub("\t", " ", CottonWorkersMouseModel.high.95.string)
CottonWorkersMouseModel.high.95.string<-tolower(CottonWorkersMouseModel.high.95.string)
CottonWorkersMouseModel.high.95.table<-table(strsplit(CottonWorkersMouseModel.high.95.string, " "))
length(CottonWorkersMouseModel.high.95)
head(sort(CottonWorkersMouseModel.high.95.table, decreasing = TRUE),60)


# now for 98

CottonWorkersMouseModel.high.98<-colnames(ReNorm.matrix)[ReNorm.matrix["CottonWorkersMouseModel",] > 0.98]
CottonWorkersMouseModel.high.98.desc<-gsm.char[match(CottonWorkersMouseModel.high.98, gsm.char$gsm),]
# how many of these contain the word lung?
length(union(grep("lung", CottonWorkersMouseModel.high.98.desc$characteristics_ch1), grep("Lung", CottonWorkersMouseModel.high.98.desc$characteristics_ch1)))


CottonWorkersMouseModel.high.98.string<-""
for (i in 1:length(unlist(CottonWorkersMouseModel.high.98.desc$characteristics_ch1))){
	CottonWorkersMouseModel.high.98.string<-paste(CottonWorkersMouseModel.high.98.string, unlist(CottonWorkersMouseModel.high.98.desc$characteristics_ch1)[i], sep = " ")
	}
CottonWorkersMouseModel.high.98.string<-gsub("[[:punct:]]", " ", CottonWorkersMouseModel.high.98.string)
CottonWorkersMouseModel.high.98.string<-gsub("\t", " ", CottonWorkersMouseModel.high.98.string)
CottonWorkersMouseModel.high.98.string<-tolower(CottonWorkersMouseModel.high.98.string)
CottonWorkersMouseModel.high.98.table<-table(strsplit(CottonWorkersMouseModel.high.98.string, " "))
length(CottonWorkersMouseModel.high.98)
head(sort(CottonWorkersMouseModel.high.98.table, decreasing = TRUE),60)

write.table(CottonWorkersMouseModel.high.95.desc, file = "/Users/GabrielAltschuler/Documents/Projects/CottonWorkers/GEOarrays_cottonHigh_0.95.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(CottonWorkersMouseModel.high.98.desc, file = "/Users/GabrielAltschuler/Documents/Projects/CottonWorkers/GEOarrays_cottonHigh_0.98.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(CottonWorkersMouseModel.high.desc, file = "/Users/GabrielAltschuler/Documents/Projects/CottonWorkers/GEOarrays_cottonHigh_1.txt", quote = FALSE, row.names = FALSE, sep = "\t")



save.image("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Data/analysisOfAdditionalDataSets.RData")


