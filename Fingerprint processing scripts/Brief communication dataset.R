# Analysis of stem cell samples using 7_10_10 build of pathway fingerprint, with 0.95, 0.05 cutoffs

setwd("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Data")
load("brief communication dataset.RData")

load("ReNorm.matrix.ternary7_10_10.RData")
database<-read.delim("GEO_ES_iPS_non_stem.txt")

# look for arrays that are not represented
database[!(database$GSM %in% colnames(ReNorm.matrix.ternary)),]
# seems to be all the MOE430B and HG-U133B - these chips were not properly represented pathway-wise so were not built into this representation. Still have 495 arrays left

database.valid<-database[(database$GSM %in% colnames(ReNorm.matrix.ternary)),]

fingerprint.database<-ReNorm.matrix.ternary[,match(database.valid$GSM, colnames(ReNorm.matrix.ternary))]

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

# remove pathways with no information content

fingerprint.database.1<-fingerprint.database[!(rowSums(fingerprint.database == 0) == length(fingerprint.database[1,])),]
fingerprint.database.1<-fingerprint.database.1[!(rowSums(fingerprint.database.1 == -1) == length(fingerprint.database.1[1,])),]
fingerprint.database.1<-fingerprint.database.1[!(rowSums(fingerprint.database.1 == 1) == length(fingerprint.database.1[1,])),]

# now only 442 pathways
# heatmap
heatmap(fingerprint.database.1, scale = "none", labCol = database.valid$CellType, margins = c(10,0), cexRow = 0.1, cexCol = 1)

database.pca<-prcomp(t(fingerprint.database))

pch = "."
cex = 4

database.celltype.colormap<-c(1:length(levels(factor(database.valid$SimpleCellType))))
names(database.celltype.colormap)<-levels(factor(database.valid$SimpleCellType))
database.celltype.colors<-database.celltype.colormap[database.valid$SimpleCellType]
# fix 


database.species.colormap<-c(1:length(levels(factor(database.valid$Species))))
names(database.species.colormap)<-levels(factor(database.valid$Species))
database.species.colors<-database.species.colormap[database.valid$Species]

plot(database.pca$x, pch = database.species.colors, cex = 1, col = database.celltype.colors)
plot(database.pca$x[,2], database.pca$x[,8], pch = database.species.colors, cex = 1, col = database.celltype.colors)
legend(-4,3, names(database.celltype.colormap), text.col = database.celltype.colormap, cex = 1, bty = "n")
legend(-3,3, names(database.species.colormap), pch = database.species.colormap, cex = 1, bty = "n")

# try just for human data
database.human<-database.valid[grep("Human", database.valid$Species),]

# remove outlier arrays
outliers<-c("GSM423940", "GSM423941", "GSM423942", "GSM423943", "GSM423944", "GSM423945")
database.human<-database.human[!(database.human$GSM %in% outliers),]

fingerprint.database.human<-ReNorm.matrix.ternary[,match(database.human$GSM, colnames(ReNorm.matrix.ternary))]
fingerprint.database.human.1<-fingerprint.database.human[!(rowSums(fingerprint.database.human == 0) == length(fingerprint.database.human[1,])),]
fingerprint.database.human.1<-fingerprint.database.human.1[!(rowSums(fingerprint.database.human.1 == -1) == length(fingerprint.database.human.1[1,])),]
fingerprint.database.human.1<-fingerprint.database.human.1[!(rowSums(fingerprint.database.human.1 == 1) == length(fingerprint.database.human.1[1,])),]
heatmap(fingerprint.database.human.1, scale = "none", labCol = database.human$CellType, margins = c(10,0), cexRow = 0.1, cexCol = 1)

database.human.pca<-prcomp(t(fingerprint.database.human))
database.human.celltype.colormap<-c(1:length(levels(factor(database.human$SimpleCellType))))
names(database.human.celltype.colormap)<-levels(factor(database.human$SimpleCellType))
database.human.celltype.colors<-database.human.celltype.colormap[database.human$SimpleCellType]

database.human.species.colormap<-c(1:length(levels(factor(database.human$Species))))
names(database.human.species.colormap)<-levels(factor(database.human$Species))
database.human.species.colors<-database.human.species.colormap[database.human$Species]

plot(database.human.pca$x, pch = database.human.species.colors, cex = 1, col = database.human.celltype.colors)

legend(-4,3, names(database.human.celltype.colormap), text.col = database.human.celltype.colormap, cex = 1, bty = "n")
legend(-3,3, names(database.human.species.colormap), pch = database.human.species.colormap, cex = 1, bty = "n")

text(database.human.pca$x, labels=database.human$GSM, cex = 0.5)

plot(database.human.pca$x[,1], database.human.pca$x[,4], pch = database.human.species.colors, cex = 1, col = database.human.celltype.colors)
text(database.human.pca$x[,1], database.human.pca$x[,4], labels=database.human$GSM, cex = 0.5)

# try extract pluripotency signature

source("/Users/GabrielAltschuler/Documents/Databases/gabriel functions.R")

pluripotent.human<-cbind(fingerprint.database.human.1[,grep("ES", database.human$SimpleCellType)], fingerprint.database.human.1[,grep("^iPS", database.human$SimpleCellType)])
testisiPS.human<-fingerprint.database.human.1[,grep("Testis-iPS", database.human$SimpleCellType)]
oligopotent.human<-fingerprint.database.human.1[,grep("Non-stem", database.human$SimpleCellType)]
pluripotent.human.sig<-rowMeans(pluripotent.human)
oligopotent.human.sig<-rowMeans(oligopotent.human)

pluri_oligo<-data.frame(pluripotent = pluripotent.human.sig, oligopotent = oligopotent.human.sig)
pluri_vs_oligo<-pluripotent.human.sig-oligopotent.human.sig
hist(abs(pluri_vs_oligo), nclass = 20)

heatmap(as.matrix(pluri_oligo)[abs(pluri_vs_oligo)>0,], scale = "none", cexCol = 1, margins = c(10,15))

pluri_vs_oligo_vs_testis.data<-database.human[match(c(colnames(pluripotent.human), colnames(oligopotent.human), colnames(testisiPS.human)), database.human$GSM),]
library(gplots)

heatmap(as.matrix(fingerprint.database.human.1)[abs(pluri_vs_oligo)>0.5,match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexCol = 1, labCol = pluri_vs_oligo_vs_testis.data$GSE, Colv = NA, margins = c(10,15))

# next step might be to try figure out why there are gaps in the heatmap - can these be correlated to storage conditions or somthing along these lines..?
# try cluster just pluripotents

heatmap(as.matrix(fingerprint.database.human.1)[abs(pluri_vs_oligo)>0.5,match(colnames(pluripotent.human), colnames(fingerprint.database.human.1))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$GSE, margins = c(10,0))

# seems like there are certain fingerprints that have Cell cycle, DNA replication and Tight junction not enriched

selectPathways<-c(grep("Cell cycle", rownames(fingerprint.database.human.1)),
					grep("DNA replication", rownames(fingerprint.database.human.1)),
					grep("Tight junction", rownames(fingerprint.database.human.1))
					)

pluripotent.subset<-colSums(fingerprint.database.human.1[selectPathways,match(colnames(pluripotent.human), colnames(fingerprint.database.human.1))])
pluripotent.subset<-pluripotent.subset[pluripotent.subset == 0]

pluri_vs_oligo_vs_testis.data[match(names(pluripotent.subset), pluri_vs_oligo_vs_testis.data$GSM),-4]

# interestingly there seems to be few false positives in the stem cell pathways (i.e. few of the fibroblast or non-stem cells demonstrate any representation in the stem cell pathways) this implies that there may be a benefit to reducing the thresholds to 10% - currently at 5% - or even going for a 2,1,0,-1,-2 representation.
# On balance, would rather have a representation where you had full coverage of 1s in the stem cell types, allowing a few 1s in the non-stem cells as you would be controlling for false positives through the number of pathways rather than single pathways.

# Processed the lower cutoffs on the server
setwd("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Data")
load("ReNorm.matrix.10ternary_7_10_10.RData")
load("ReNorm.matrix.15ternary_7_10_10.RData")
load("ReNorm.matrix.20ternary_7_10_10.RData")
load("ReNorm.matrix.25ternary_7_10_10.RData")

fingerprint.database.human.10<-ReNorm.matrix.ternary.10[,match(database.human$GSM, colnames(ReNorm.matrix.ternary.10))]
fingerprint.database.human.15<-ReNorm.matrix.ternary.15[,match(database.human$GSM, colnames(ReNorm.matrix.ternary.15))]
fingerprint.database.human.20<-ReNorm.matrix.ternary.20[,match(database.human$GSM, colnames(ReNorm.matrix.ternary.20))]
fingerprint.database.human.25<-ReNorm.matrix.ternary.25[,match(database.human$GSM, colnames(ReNorm.matrix.ternary.25))]

temp<-heatmap(fingerprint.database.human[names(pluri_vs_oligo[abs(pluri_vs_oligo)>0.5]),match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15))

heatmap(fingerprint.database.human.10[names(pluri_vs_oligo[abs(pluri_vs_oligo)>0.5]),match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15))

heatmap(fingerprint.database.human.15[names(pluri_vs_oligo[abs(pluri_vs_oligo)>0.5]),match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15))

heatmap(fingerprint.database.human.20[names(pluri_vs_oligo[abs(pluri_vs_oligo)>0.5]),match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15))

heatmap(fingerprint.database.human.25[names(pluri_vs_oligo[abs(pluri_vs_oligo)>0.5]),match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15))


# PCA

database.human.pca.10<-prcomp(t(fingerprint.database.human.10))
database.human.pca.15<-prcomp(t(fingerprint.database.human.15))
database.human.pca.20<-prcomp(t(fingerprint.database.human.20))
plot(database.human.pca$x[,1], database.human.pca$x[,2], pch = database.human.species.colors, cex = 1, col = database.human.celltype.colors)
plot(database.human.pca.10$x[,1], database.human.pca.10$x[,2], pch = database.human.species.colors, cex = 1, col = database.human.celltype.colors)
plot(database.human.pca.15$x[,1], database.human.pca.15$x[,2], pch = database.human.species.colors, cex = 1, col = database.human.celltype.colors)
plot(database.human.pca.20$x[,1], database.human.pca.20$x[,2], pch = database.human.species.colors, cex = 1, col = database.human.celltype.colors)

legend(-4,3, names(database.human.celltype.colormap), text.col = database.human.celltype.colormap, cex = 1, bty = "n")
legend(-3,3, names(database.human.species.colormap), pch = database.human.species.colormap, cex = 1, bty = "n")


fingerprint.10.database<-ReNorm.matrix.ternary.10[,match(database.valid$GSM, colnames(ReNorm.matrix.ternary.10))]
fingerprint.15.database<-ReNorm.matrix.ternary.15[,match(database.valid$GSM, colnames(ReNorm.matrix.ternary.15))]
fingerprint.20.database<-ReNorm.matrix.ternary.20[,match(database.valid$GSM, colnames(ReNorm.matrix.ternary.20))]
fingerprint.25.database<-ReNorm.matrix.ternary.25[,match(database.valid$GSM, colnames(ReNorm.matrix.ternary.25))]

database.pca<-prcomp(t(fingerprint.database))
database.10.pca<-prcomp(t(fingerprint.10.database))
database.15.pca<-prcomp(t(fingerprint.15.database))
database.20.pca<-prcomp(t(fingerprint.20.database))
database.25.pca<-prcomp(t(fingerprint.25.database))

plot(database.pca$x, pch = database.species.colors, cex = 1, col = database.celltype.colors)
plot(database.10.pca$x, pch = database.species.colors, cex = 1, col = database.celltype.colors)
plot(database.15.pca$x, pch = database.species.colors, cex = 1, col = database.celltype.colors)
plot(database.20.pca$x, pch = database.species.colors, cex = 1, col = database.celltype.colors)
plot(database.25.pca$x, pch = database.species.colors, cex = 1, col = database.celltype.colors)

# re-assess signiture at higher threshold

pluripotent.human.10<-cbind(fingerprint.database.human.10[,grep("ES", database.human$SimpleCellType)], fingerprint.database.human.10[,grep("^iPS", database.human$SimpleCellType)])

oligopotent.human.10<-fingerprint.database.human.10[,grep("Non-stem", database.human$SimpleCellType)]
pluripotent.human.10.sig<-rowMeans(pluripotent.human.10)
oligopotent.human.10.sig<-rowMeans(oligopotent.human.10)

pluri_oligo.10<-data.frame(pluripotent = pluripotent.human.10.sig, oligopotent = oligopotent.human.10.sig)
pluri_vs_oligo.10<-pluripotent.human.10.sig-oligopotent.human.10.sig

pluripotent.human.15<-cbind(fingerprint.database.human.15[,grep("ES", database.human$SimpleCellType)], fingerprint.database.human.15[,grep("^iPS", database.human$SimpleCellType)])

oligopotent.human.15<-fingerprint.database.human.15[,grep("Non-stem", database.human$SimpleCellType)]
pluripotent.human.15.sig<-rowMeans(pluripotent.human.15)
oligopotent.human.15.sig<-rowMeans(oligopotent.human.15)

pluri_oligo.15<-data.frame(pluripotent = pluripotent.human.15.sig, oligopotent = oligopotent.human.15.sig)
pluri_vs_oligo.15<-pluripotent.human.15.sig-oligopotent.human.15.sig



pluripotent.human.20<-cbind(fingerprint.database.human.20[,grep("ES", database.human$SimpleCellType)], fingerprint.database.human.20[,grep("^iPS", database.human$SimpleCellType)])

oligopotent.human.20<-fingerprint.database.human.20[,grep("Non-stem", database.human$SimpleCellType)]
pluripotent.human.20.sig<-rowMeans(pluripotent.human.20)
oligopotent.human.20.sig<-rowMeans(oligopotent.human.20)

pluri_oligo.20<-data.frame(pluripotent = pluripotent.human.20.sig, oligopotent = oligopotent.human.20.sig)
pluri_vs_oligo.20<-pluripotent.human.20.sig-oligopotent.human.20.sig

pluripotent.human.25<-cbind(fingerprint.database.human.25[,grep("ES", database.human$SimpleCellType)], fingerprint.database.human.25[,grep("^iPS", database.human$SimpleCellType)])

oligopotent.human.25<-fingerprint.database.human.25[,grep("Non-stem", database.human$SimpleCellType)]
pluripotent.human.25.sig<-rowMeans(pluripotent.human.25)
oligopotent.human.25.sig<-rowMeans(oligopotent.human.25)

pluri_oligo.25<-data.frame(pluripotent = pluripotent.human.25.sig, oligopotent = oligopotent.human.25.sig)
pluri_vs_oligo.25<-pluripotent.human.25.sig-oligopotent.human.25.sig




# compare heatmaps of the top plutipotency pathways derived from 15% threshold data

heatmap(fingerprint.database.human[names(pluri_vs_oligo.15[abs(pluri_vs_oligo.15)>0.5]),match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15))

heatmap(fingerprint.database.human.15[names(pluri_vs_oligo.15[abs(pluri_vs_oligo.15)>0.5]),match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15))

# seems that the 15% data might work well - should probably think of a way of defining this properly - perhaps an average distance between stem and non-stem cell types # does this need to be divided by something..? perhaps the abs total for each - not sure this is correct as highest ratio would be for 
# decided on the length of the difference vector divided by the leng

sqrt(sum(pluri_vs_oligo^2))/(sqrt(sum((pluri_oligo^2)[,1]))-sqrt(sum((pluri_oligo^2)[,2])))
sqrt(sum(pluri_vs_oligo.10^2))/(sqrt(sum((pluri_oligo.10^2)[,1]))-sqrt(sum((pluri_oligo.10^2)[,2])))
sqrt(sum(pluri_vs_oligo.15^2))/(sqrt(sum((pluri_oligo.15^2)[,1]))-sqrt(sum((pluri_oligo.15^2)[,2])))
sqrt(sum(pluri_vs_oligo.20^2))/(sqrt(sum((pluri_oligo.20^2)[,1]))-sqrt(sum((pluri_oligo.20^2)[,2])))
sqrt(sum(pluri_vs_oligo.25^2))/(sqrt(sum((pluri_oligo.25^2)[,1]))-sqrt(sum((pluri_oligo.25^2)[,2])))

sqrt(sum(pluri_vs_oligo^2))
sqrt(sum(pluri_vs_oligo.10^2))
sqrt(sum(pluri_vs_oligo.15^2))
sqrt(sum(pluri_vs_oligo.20^2))
sqrt(sum(pluri_vs_oligo.25^2))

library(gplots)
heatmap.2(fingerprint.database.human.15[names(pluri_vs_oligo.15[abs(pluri_vs_oligo.15)>0.75]),match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexRow = 0.5, cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15), trace = "none", dendrogram = "row", col = colorRampPalette(c("blue","grey", "red")))

# cross-species comparison

pluripotent.10<-cbind(fingerprint.10.database[,grep("ES", database.valid$SimpleCellType)], fingerprint.10.database[,grep("^iPS", database.valid$SimpleCellType)])

testisiPS.10<-fingerprint.10.database[,grep("Testis-iPS", database.valid$SimpleCellType)]
testisiPS.15<-fingerprint.15.database[,grep("Testis-iPS", database.valid$SimpleCellType)]
testisiPS.25<-fingerprint.15.database[,grep("Testis-iPS", database.valid$SimpleCellType)]

pluri_vs_oligo_vs_testis.xspecies.data<-database.valid[match(c(colnames(pluripotent.10), colnames(oligopotent.10), colnames(testisiPS.10)), database.valid$GSM),]



oligopotent.10<-fingerprint.10.database[,grep("Non-stem", database.valid$SimpleCellType)]
pluripotent.10.sig<-rowMeans(pluripotent.10)
oligopotent.10.sig<-rowMeans(oligopotent.10)

pluri_oligo.xspecies.10<-data.frame(pluripotent = pluripotent.10.sig, oligopotent = oligopotent.10.sig)
pluri_vs_oligo.xspecies.10<-pluripotent.10.sig-oligopotent.10.sig

heatmap.2(fingerprint.10.database[names(pluri_vs_oligo.xspecies.10[abs(pluri_vs_oligo.xspecies.10)>0.75]),match(pluri_vs_oligo_vs_testis.xspecies.data$GSM, colnames(fingerprint.10.database))], scale = "none", cexCol = 0.5, labCol = pluri_vs_oligo_vs_testis.xspecies.data$Species, Colv = NA, margins = c(10,15), trace = "none", dendrogram = "row", col = colorRampPalette(c("blue","grey", "red")))

# remove pathways that are 
final.pathways<-names(pluri_vs_oligo.15[abs(pluri_vs_oligo.15)>0.75])
final.pathways<-final.pathways[!(grepl("\\{", final.pathways))]
final.pathways<-final.pathways[!(grepl("Blakely", final.pathways))]
final.pathways.names<-final.pathways
# update names

# update netpaths TR gene sets
final.pathways.names[grep("DOWN", final.pathways.names)]<-paste("Down reg. targets of", final.pathways.names[grep("DOWN", final.pathways.names)], sep = " ")
final.pathways.names[grep("UP", final.pathways.names)]<-paste("Up reg. targets of", final.pathways.names[grep("UP", final.pathways.names)], sep = " ")
final.pathways.names[grep("DOWN", final.pathways.names)]<-gsub("DOWN", "sig. pathway", final.pathways.names)[grep("DOWN", final.pathways.names)]
final.pathways.names[grep("UP", final.pathways.names)]<-gsub("UP", "sig. pathway", final.pathways.names)[grep("UP", final.pathways.names)]
# Correct KEGG names to remove number string
final.pathways.names[grep("^0", final.pathways.names)]<-substr(final.pathways.names, 7, 200)[grep("^0", final.pathways.names)]

# capitalize first word
startCap<-function(x){
	paste(toupper(substring(x, 1,1)), substring(x, 2),sep="", collapse=" ")
	}

final.pathways.names<-unlist(lapply(final.pathways.names, startCap))

# difficult to decapitalise certain words and not others as get things like dna etc.
final.pathways.names<-gsub("Receptor", "receptor", final.pathways.names)
final.pathways.names<-gsub("Pathway", "pathway", final.pathways.names)
final.pathways.names<-gsub("Cascade", "cascade", final.pathways.names)
final.pathways.names<-gsub("Signaling", "signaling", final.pathways.names)
final.pathways.names<-gsub("Selenoproteins", "selenoproteins", final.pathways.names)
final.pathways.names<-gsub("Synthesis", "synthesis", final.pathways.names)
final.pathways.names<-gsub("Transcription", "transcription", final.pathways.names)
final.pathways.names<-gsub("Initiation", "initiation", final.pathways.names)
final.pathways.names<-gsub("Biosynthesis", "biosynthesis", final.pathways.names)
final.pathways.names<-gsub("Proteins", "proteins", final.pathways.names)
final.pathways.names<-gsub("Degradation", "degradation", final.pathways.names)
final.pathways.names<-gsub("Ribosomal", "ribosomal", final.pathways.names)
final.pathways.names<-gsub("Replication", "replication", final.pathways.names)
final.pathways.names<-gsub("Acid", "acid", final.pathways.names)
final.pathways.names<-gsub("Id", "ID", final.pathways.names)
final.pathways.names<-unlist(lapply(final.pathways.names, startCap))
final.pathways.names<-gsub("P53", "p53", final.pathways.names)

heatmap.2(fingerprint.database.human.15[final.pathways,match(pluri_vs_oligo_vs_testis.data$GSM, colnames(fingerprint.database.human.1))], scale = "none", cexRow = 0.5, cexCol = 0.5, labRow = final.pathways.names, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15), trace = "none", dendrogram = "row", col = colorRampPalette(c("blue","grey", "red")))

# remove outlier GSM248213
final.GSM<-pluri_vs_oligo_vs_testis.data$GSM[!(grepl("GSM248213", pluri_vs_oligo_vs_testis.data$GSM))]
final.CellType<-pluri_vs_oligo_vs_testis.data$SimpleCellType[!(grepl("GSM248213", pluri_vs_oligo_vs_testis.data$GSM))]

heatmap.2(fingerprint.database.human.15[final.pathways,match(final.GSM, colnames(fingerprint.database.human.1))], scale = "none", cexRow = 1, cexCol = 0.5, labRow = final.pathways.names, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,20), trace = "none", dendrogram = "row", col = colorRampPalette(c("blue","grey", "red")))

# full fingerprint
heatmap.2(fingerprint.database.human[,match(final.GSM, colnames(fingerprint.database.human.1))], scale = "none", cexRow = 0.5, cexCol = 0.5, labRow = final.pathways.names, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15), trace = "none", dendrogram = "row", col = colorRampPalette(c("blue","grey", "red")))
heatmap.2(fingerprint.database.human.10[,match(final.GSM, colnames(fingerprint.database.human.1))], scale = "none", cexRow = 0.5, cexCol = 0.5, labRow = final.pathways.names, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15), trace = "none", dendrogram = "row", col = colorRampPalette(c("blue","grey", "red")))
heatmap.2(fingerprint.database.human.15[,match(final.GSM, colnames(fingerprint.database.human.1))], scale = "none", cexRow = 0.5, cexCol = 0.5, labRow = final.pathways.names, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15), trace = "none", dendrogram = "row", col = colorRampPalette(c("blue","grey", "red")))
heatmap.2(fingerprint.database.human.20[,match(final.GSM, colnames(fingerprint.database.human.1))], scale = "none", cexRow = 0.5, cexCol = 0.5, labRow = final.pathways.names, labCol = pluri_vs_oligo_vs_testis.data$SimpleCellType, Colv = NA, margins = c(10,15), trace = "none", dendrogram = "row", col = colorRampPalette(c("blue","grey", "red")))



# perhaps might be nice to plot an overlap matrix of the stem cell pathways - can use previously calculated network from cytoscape - just need to add an attribute file for stem levels

pathways<-rownames(fingerprint.database)
pathway.names<-rownames(fingerprint.database)

pathway.names[grep("DOWN", pathway.names)]<-paste("Down reg. targets of", pathway.names[grep("DOWN", pathway.names)], sep = " ")
pathway.names[grep("UP", pathway.names)]<-paste("Up reg. targets of", pathway.names[grep("UP", pathway.names)], sep = " ")
pathway.names[grep("DOWN", pathway.names)]<-gsub("DOWN", "sig. pathway", pathway.names)[grep("DOWN", pathway.names)]
pathway.names[grep("UP", pathway.names)]<-gsub("UP", "sig. pathway", pathway.names)[grep("UP", pathway.names)]
# Correct KEGG names to remove number string
pathway.names[grep("^0", pathway.names)]<-substr(pathway.names, 7, 200)[grep("^0", pathway.names)]

# capitalize first word
startCap<-function(x){
	paste(toupper(substring(x, 1,1)), substring(x, 2),sep="", collapse=" ")
	}

pathway.names<-unlist(lapply(pathway.names, startCap))

# difficult to decapitalise certain words and not others as get things like dna etc.
pathway.names<-gsub("Receptor", "receptor", pathway.names)
pathway.names<-gsub("Pathway", "pathway", pathway.names)
pathway.names<-gsub("Cascade", "cascade", pathway.names)
pathway.names<-gsub("Signaling", "signaling", pathway.names)
pathway.names<-gsub("Selenoproteins", "selenoproteins", pathway.names)
pathway.names<-gsub("Synthesis", "synthesis", pathway.names)
pathway.names<-gsub("Transcription", "transcription", pathway.names)
pathway.names<-gsub("Initiation", "initiation", pathway.names)
pathway.names<-gsub("Biosynthesis", "biosynthesis", pathway.names)
pathway.names<-gsub("Proteins", "proteins", pathway.names)
pathway.names<-gsub("Degradation", "degradation", pathway.names)
pathway.names<-gsub("Ribosomal", "ribosomal", pathway.names)
pathway.names<-gsub("Replication", "replication", pathway.names)
pathway.names<-gsub("Acid", "acid", pathway.names)
pathway.names<-gsub("Id", "ID", pathway.names)
pathway.names<-unlist(lapply(pathway.names, startCap))
pathway.names<-gsub("P53", "p53", pathway.names)


pathways.frame<-data.frame(Name = pathways, NameUpdate = pathway.names, pluri = 100*pluri_vs_oligo.15)
write.table(pathways.frame, file = "pluri_vs_oligo.15.txt", quote = FALSE, row.names = FALSE, sep = "\t")

# An interesting aspect would be to examine the PCA of the entire GEO space to see where the cells lie...
# this doesn't provide anything particularly interesting

full.pca<-prcomp(t(ReNorm.matrix.ternary))
full.pca.human<-full.pca$x[database.human$GSM,]
plot(full.pca.human[,1], full.pca.human[,2], pch = database.human.species.colors, cex = 1, col = database.human.celltype.colors)

# plot simple pca for figure
database.15.pca.figure<-prcomp(t(fingerprint.15.database[,match(final.GSM, colnames(fingerprint.15.database))]))

database.15.pca.figure.colormap<-c("magenta","magenta","brown","green")
names(database.15.pca.figure.colormap)<-levels(factor(final.CellType))
database.15.pca.figure.colors<-database.15.pca.figure.colormap[final.CellType]

plot(database.15.pca.figure$x[,1], database.15.pca.figure$x[,2], pch = 16, cex = 1, col = database.15.pca.figure.colors, ylim = c(-7.5,7.5), xlim = c(-12,8))

# could use the pluripotency fingerprint to 'fish out' similar GEO samples

# First define pluripotency fingerprint
pluri.fingerprint<-rep(0,length(names(pluri_vs_oligo.15)))
names(pluri.fingerprint)<-names(pluri_vs_oligo.15)
pluri.fingerprint[pluri_vs_oligo.15>0.75]<-1
pluri.fingerprint[pluri_vs_oligo.15<(-0.75)]<-(-1)
# pluri.fingerprint<-as.data.frame(pluri.fingerprint)

# now subtract this fingerprint from the full database and subset for the stem pathways
pluri.sub<-ReNorm.matrix.ternary.15-pluri.fingerprint
pluri.sub<-pluri.sub[(abs(pluri_vs_oligo.15)>0.75),]

# subtract pluri.fingerprint from full GEO corpus and see which colsum is the minimum
pluri.diff<-colSums(abs(pluri.sub))
head(pluri.diff[order(pluri.diff)])

# histogram of differences
hist(pluri.diff, nclass = 100)
hist(pluri.diff, nclass = 500, xlim = c(0,20), ylim = c(0,200))

# now need to think about how to benchmark these samples

# take sample of half of the pluripotent cells. Build signature, don't bother with subtraction of non-stem cell set as should already be normalized by virtue of the cutoff. Then take the abs > 0.75 to build the stem cell signature. Use this to 1) build the GEO set (slow) and 2) determine the distances of the remaining half, in data set have 216 pluripotent and 56 oligopotent samples.

thresholded.fingerprints<-fingerprint.database.human.15
signature.threshold<-0.75
thresholded.fingerprints.full<-cd133.GEO.ternary.15

# define pluripotent (216) and oligopotent (57) samples
pluripotents<-as.character(database.human[c(grep("ES", database.human$SimpleCellType) ,grep("^iPS", database.human$SimpleCellType)),3])
oligopotents<-as.character(database.human[grep("Non-stem", database.human$SimpleCellType),3])
testis<-as.character(database.human[grep("Testis-iPS", database.human$SimpleCellType),3])
# remove outlier
oligopotents<-oligopotents[!(grepl("GSM248213", oligopotents))]
pluripotent.distance.record<-rep(0, length(thresholded.fingerprints.full[1,]))
names(pluripotent.distance.record)<-colnames(thresholded.fingerprints.full)

# start loop here
for (i in 1:100){
	pluripotent.test.sample<-sample(pluripotents, round(0.99*length(pluripotents),0))
	pluripotent.control.sample<-pluripotents[!(pluripotents %in% pluripotent.test.sample)]

	pluripotent.test.fingerprints<-thresholded.fingerprints[,match(pluripotent.test.sample, colnames(thresholded.fingerprints))]
	pluripotent.control.fingerprints<-thresholded.fingerprints[,match(pluripotent.control.sample, colnames(thresholded.fingerprints))]

	pluripotent.test.signature<-rowMeans(pluripotent.test.fingerprints)
	pluripotent.test.signature.fingerprint<-rep(0,length(names(pluripotent.test.signature)))
	names(pluripotent.test.signature.fingerprint)<-names(pluripotent.test.signature)
	pluripotent.test.signature.fingerprint[pluripotent.test.signature>signature.threshold]<-1
	pluripotent.test.signature.fingerprint[pluripotent.test.signature<(-signature.threshold)]<-(-1)

	# build full distribution
	pluripotent.difference<-thresholded.fingerprints.full[(abs(pluripotent.test.signature)>signature.threshold),] - pluripotent.test.signature.fingerprint[(abs(pluripotent.test.signature)>signature.threshold)]
	pluripotent.distance<-colSums(abs(pluripotent.difference))
	pluripotent.distance.record<-pluripotent.distance+pluripotent.distance.record
	}
pluripotent.distance.record<-pluripotent.distance.record/100
pdf("100xfull_sample signature_15_0.9.pdf")
hist(pluripotent.distance, xlim = c(0,200), breaks = seq(0,200,1), main = "Full data corpus", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[pluripotent.test.sample], xlim = c(0,200), breaks = seq(0,200,1), main = "ES/iPS test sample used to build signature", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[pluripotent.control.sample], xlim = c(0,200), breaks = seq(0,200,1), main = "ES/iPS control, not used to build signature", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[pluripotents], xlim = c(0,200), breaks = seq(0,200,1), main = "ES/iPS sample", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[oligopotents], xlim = c(0,200), breaks = seq(0,200,1), main = "Curated set of non-stem cells", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[pluripotent.test.sample], xlim = c(0,200), breaks = seq(0,200,1), main = "ES/iPS test sample used to build signature, testis-iPS in red", xlab = "distance from pluripotency signature", mar = c(0,20))
hist(pluripotent.distance.record[testis], xlim = c(0,200), breaks = seq(0,200,1), col ="red", add = TRUE)
par(mfcol = c(2,1))
hist(pluripotent.distance, xlim = c(0,200), breaks = seq(0,200,1), main = "Full data corpus", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[pluripotents], xlim = c(0,200), breaks = seq(0,200,1), main = "ES/iPS sample", xlab = "distance from pluripotency signature", col = "blue")
hist(pluripotent.distance.record[testis], xlim = c(0,200), breaks = seq(0,200,1), col ="red", add = TRUE)
dev.off()

# Use median or mean rank to asses the pluripotency distance measure (averaged over 100 samples)

# run with cd133 data
load("/Users/GabrielAltschuler/Documents/Projects/Melanoma/cd133.GEO.ternary.15.RData")
cd133.pos<-colnames(cd133.GEO.ternary.15)[grep("POS", colnames(cd133.GEO.ternary.15))]
cd133.neg<-colnames(cd133.GEO.ternary.15)[grep("NEG", colnames(cd133.GEO.ternary.15))]
hist(pluripotent.distance.record[cd133.pos], xlim = c(0,200), breaks = seq(0,200,1), main = "CD133 posative", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[cd133.neg], xlim = c(0,200), breaks = seq(0,200,1), main = "CD133 negative", xlab = "distance from pluripotency signature")



# compare to ES in GPL571
GPL571.pluri<-c("GSM339011", "GSM339012", "GSM339425", "GSM339426", "GSM339427", "GSM424314", "GSM424315", 
"GSM424316","GSM424317","GSM424318","GSM424319")

GPL571.pluri1<-c("GSM339011","GSM339425", "GSM339427")

hist(pluripotent.distance.record[GPL571.pluri], xlim = c(0,200), breaks = seq(0,200,1), main = "differentiated ES", xlab = "distance from pluripotency signature")
hoxa9<-c("GSM344801", "GSM344802", "GSM344803", "GSM344804", "GSM344805", "GSM344806") Control<-c("GSM344807", "GSM344808", "GSM344809", "GSM344810", "GSM344811", "GSM344812")

pdf("/Users/GabrielAltschuler/Documents/Projects/Melanoma/pluripotent distance.pdf")
hist(pluripotent.distance, xlim = c(0,200), breaks = seq(0,200,1), main = "Full data corpus", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[pluripotents], xlim = c(0,200), breaks = seq(0,200,1), main = "ES/iPS sample", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[oligopotents], xlim = c(0,200), breaks = seq(0,200,1), main = "Curated set of non-stem cells", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[cd133.pos], xlim = c(0,200), breaks = seq(0,200,1), main = "CD133 positive", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[cd133.neg], xlim = c(0,200), breaks = seq(0,200,1), main = "CD133 negative", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[GPL571.pluri], xlim = c(0,200), breaks = seq(0,200,1), main = "differentiated ES", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[hoxa9], xlim = c(0,200), breaks = seq(0,200,1), main = "Armstrong 2009 - HoxA9 MLL", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[Control], xlim = c(0,200), breaks = seq(0,200,1), main = "Armstrong 2009 - MLL control", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[blood[grep("HSC", names(blood))]], xlim = c(0,200), breaks = seq(0,200,1), main = "HSC", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[blood[grep("MEP", names(blood))]], xlim = c(0,200), breaks = seq(0,200,1), main = "MEP", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[blood[grep("GMP", names(blood))]], xlim = c(0,200), breaks = seq(0,200,1), main = "GMP", xlab = "distance from pluripotency signature")
hist(pluripotent.distance.record[blood[grep("CMP", names(blood))]], xlim = c(0,200), breaks = seq(0,200,1), main = "CMP", xlab = "distance from pluripotency signature")
dev.off()

# look for closest fingerprint

neg = rowMeans(cd133.ternary.15[,grep("NEG", colnames(cd133.ternary.15))])
pos = rowMeans(cd133.ternary.15[,grep("POS", colnames(cd133.ternary.15))])

cd133.neg.order<-cd133.neg[c(3,6,1,4,5,2,7)]
cd133.pos.order<-cd133.pos[c(3,6,1,4,5,2,7)]

cd133.pathway.frame<-cd133.ternary.15[abs(neg-pos)>0.2,c(cd133.neg.order, cd133.pos.order)]
save(cd133.pathway.frame, file = "/Users/GabrielAltschuler/Documents/Projects/Melanoma/cd133.pathway.frame.RData")
pdf("/Users/GabrielAltschuler/Documents/Projects/Melanoma/cd133Heatmap.pdf")
heatmap(cd133.pathway.frame, scale = "none", Colv = NA, margins = c(7,10), cexRow = 0.5, cexCol = 0.5)
dev.off()

neg.test.signature.fingerprint<-rep(0,length(names(neg)))
names(neg.test.signature.fingerprint)<-names(neg)
neg.test.signature.fingerprint[neg>signature.threshold]<-1
neg.test.signature.fingerprint[neg<(-signature.threshold)]<-(-1)

neg.difference<-thresholded.fingerprints.full[(abs(neg)>signature.threshold),] - neg.test.signature.fingerprint[(abs(neg)>signature.threshold)]
neg.distance<-colSums(abs(neg.difference))
hist(neg.distance, xlim = c(0,100), breaks = seq(0,200,1), main = "Full data corpus", xlab = "distance from CD133 Negative signature")
hist(neg.distance, xlim = c(0,100), breaks = seq(0,200,1), main = "Full data corpus", xlab = "distance from CD133 Negative signature")
hist(neg.distance[oligopotents], xlim = c(0,100), breaks = seq(0,200,1), main = "Non-stem", xlab = "distance from CD133 Negative signature")

pos.test.signature.fingerprint<-rep(0,length(names(pos)))
names(pos.test.signature.fingerprint)<-names(pos)
pos.test.signature.fingerprint[pos>signature.threshold]<-1
pos.test.signature.fingerprint[pos<(-signature.threshold)]<-(-1)

pos.test.signature.fingerprint.subset<-pos.test.signature.fingerprint[abs(pos.test.signature.fingerprint)>0]

temp<-as.data.frame(pos.test.signature.fingerprint.subset[order(pos.test.signature.fingerprint.subset)])
colnames(temp)<-"score"
pos.difference<-thresholded.fingerprints.full[(abs(pos)>signature.threshold),] - pos.test.signature.fingerprint[(abs(pos)>signature.threshold)]
pos.distance<-colSums(abs(pos.difference))
pdf("/Users/GabrielAltschuler/Documents/Projects/Melanoma/CD133 positive distance.pdf")
hist(pos.distance, xlim = c(0,100), breaks = seq(0,200,1), main = "Full data corpus", xlab = "distance from CD133 positive signature")
hist(pos.distance[pluripotents], xlim = c(0,100), breaks = seq(0,200,1), main = "ES/iPS", xlab = "distance from CD133 positive signature")
hist(pos.distance[oligopotents], xlim = c(0,100), breaks = seq(0,200,1), main = "Non-stem", xlab = "distance from CD133 positive signature")
hist(pos.distance[cd133.pos], xlim = c(0,200), breaks = seq(0,200,1), main = "CD133 positive", xlab = "distance from CD133 positive signature")
hist(pos.distance[cd133.neg], xlim = c(0,200), breaks = seq(0,200,1), main = "CD133 negative", xlab = "distance from CD133 positive signature")
dev.off()
# Top 30 were all Human Melanoma samples, from 3 different platforms
# Expression profiling of BRAFV600E melanoma cell lines upon RAF inhibition HG-U133A_2
# In vitro study of gene expression and response to chemotherapy across 50 human-melanoma derived cell lines, HGU133Plus2
# NCI60 cancer cell line ME:SK-MEL-28 HG-U133A
# Melanoma cell line D25 HGU133Plus2
# SK-MEL-28 hg133a21
# save image upto this point (dataframe does not contain all the above data - just the subset necesary for the cd133 analysis)

# re-color for presentation

pdf("/Users/GabrielAltschuler/Documents/Projects/Melanoma/pluripotent distance - recolor.pdf")
par(mfcol = c(2,1), mar = c(0, 4, 4, 0))
hist(pluripotent.distance, xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "", margin = c(0, 4, 4, 0) )
par(mar = c(7, 4, 4, 0))
hist(pluripotent.distance.record[pluripotents], xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "distance from pluripotency signature", col = "green", margin = c(5, 4, 4, 2))
hist(pluripotent.distance.record[cd133.pos], xlim = c(0,200), breaks = seq(0,200,1), col = "red", main = "", add = TRUE)
hist(pluripotent.distance.record[cd133.neg], xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "distance from pluripotency signature", col = "blue", add = TRUE)
dev.off()


pdf("/Users/GabrielAltschuler/Documents/Projects/Melanoma/pluripotent distance - recolor - separates.pdf")
hist(pluripotent.distance.record[cd133.pos], xlim = c(0,200), breaks = seq(0,200,1), main = "CD133 positive", xlab = "distance from pluripotency signature", col = "red")
hist(pluripotent.distance.record[cd133.neg], xlim = c(0,200), breaks = seq(0,200,1), main = "CD133 negative", xlab = "distance from pluripotency signature", col = "blue")
hist(pluripotent.distance.record[GPL571.pluri], xlim = c(0,200), breaks = seq(0,200,1), main = "differentiated ES", xlab = "distance from pluripotency signature", col = "green")
dev.off()




# stats for paper

mean(pluripotent.distance.record)
max(pluripotent.distance.record[pluripotents])
sum(pluripotent.distance.record<mean(pluripotent.distance.record[pluripotents]))/sum(pluripotent.distance.record>mean(pluripotent.distance.record[pluripotents]))

sum(pluripotent.distance.record<max(pluripotent.distance.record[pluripotents]))/sum(pluripotent.distance.record>max(pluripotent.distance.record[pluripotents]))

1/sum(pluripotent.distance.record>min(pluripotent.distance.record[pluripotents]))

sum(pluripotent.distance.record<mean(pluripotent.distance.record[testis]))/sum(pluripotent.distance.record>max(pluripotent.distance.record[testis]))

# define testis signature

testis.sig = rowMeans(fingerprint.database.human.15[,testis])

testis.sig.fingerprint<-rep(0,length(names(testis.sig)))
names(testis.sig.fingerprint)<-names(testis.sig)
testis.sig.fingerprint[testis.sig>signature.threshold]<-1
testis.sig.fingerprint[testis.sig<(-signature.threshold)]<-(-1)

testis.sig.difference<-thresholded.fingerprints.full[(abs(testis.sig)>signature.threshold),] - testis.sig.fingerprint[(abs(testis.sig)>signature.threshold)]
testis.sig.distance<-colSums(abs(testis.sig.difference))


blood<-c("GSM300093", "GSM300100", "GSM300101", "GSM300102", "GSM300103", "GSM300104", "GSM300105", "GSM300106", "GSM300107", "GSM300108", "GSM300109", "GSM300111", "GSM300112", "GSM300114", "GSM300115", "GSM300116", "GSM300117", "GSM300118", "GSM300119", "GSM300120", "GSM300121", "GSM300122", "GSM300123", "GSM300124", "GSM300125", "GSM300126", "GSM300127", "GSM300128", "GSM300129", "GSM300130", "GSM300131", "GSM300132", "GSM300133", "GSM300134", "GSM300135", "GSM300136", "GSM300137", "GSM300138", "GSM300139", "GSM300140", "GSM300141", "GSM300142", "GSM300144", "GSM300146", "GSM300147", "GSM300148", "GSM300149")

names(blood)<-c("CML_1_HSC", "CML_1_GMP", "CML_1_MEP", "CML_2_HSC", "CML_2_CMP", "CML_2_GMP", "CML_2_MEP", "CML_3_HSC", "CML_3_CMP", "CML_3_GMP", "CML_3_MEP", "CML_4_HSC", "CML_4_CMP", "CML_4_GMP", "CML_4_MEP", "CML_5_HSC", "CML_5_CMP", "CML_5_GMP", "CML_5_MEP", "CML_6_HSC", "CML_6_CMP", "CML_6_GMP", "CML_6_MEP", "CML_7_HSC", "CML_7_CMP", "CML_7_GMP", "CML_7_MEP", "CTRL_1_HSC", "CTRL_1_CMP", "CTRL_1_GMP", "CTRL_1_MEP", "CTRL_2_HSC", "CTRL_2_CMP", "CTRL_2_GMP", "CTRL_2_MEP", "CTRL_3_HSC", "CTRL_3_CMP", "CTRL_3_GMP", "CTRL_3_MEP", "CTRL_4_HSC", "CTRL_4_CMP", "CTRL_4_GMP", "CTRL_4_MEP", "CTRL_5_HSC", "CTRL_5_CMP", "CTRL_5_GMP", "CTRL_5_MEP")

#######
# Analysis of Rossi mRNA derived stem cells


rossi<-read.delim("Rossi_RiPS.txt", header = FALSE)
pdf(file = "Rossi_RiPS_distance_from_pluripotency.pdf")
par(mfcol = c(2,1), mar = c(0, 4, 4, 0))

hist(pluripotent.distance, xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "", margin = c(0, 4, 4, 0) )
par(mar = c(7, 4, 4, 0))
hist(pluripotent.distance.record[pluripotents], xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "distance from pluripotency signature", col = "green", margin = c(5, 4, 4, 2))
hist(pluripotent.distance.record[as.character(rossi[,1])], xlim = c(0,200), breaks = seq(0,200,1), col = "red", main = "", add = TRUE)

par(mfcol = c(2,1), mar = c(0, 4, 4, 0))

hist(pluripotent.distance, xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "", margin = c(0, 4, 4, 0) )
par(mar = c(7, 4, 4, 0))
hist(pluripotent.distance.record[pluripotents], xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "distance from pluripotency signature", col = "green", margin = c(5, 4, 4, 2))
hist(pluripotent.distance.record[as.character(rossi[grep("dH1F", rossi[,2]),1])], xlim = c(0,200), breaks = seq(0,200,1), col = "red", main = "", add = TRUE)
hist(pluripotent.distance.record[as.character(rossi[grep("BJ", rossi[,2]),1])], xlim = c(0,200), breaks = seq(0,200,1), col = "yellow", main = "", add = TRUE)
hist(pluripotent.distance.record[as.character(rossi[grep("CF", rossi[,2]),1])], xlim = c(0,200), breaks = seq(0,200,1), col = "black", main = "", add = TRUE)
hist(pluripotent.distance.record[as.character(rossi[grep("MRC5", rossi[,2]),1])], xlim = c(0,200), breaks = seq(0,200,1), col = "blue", main = "", add = TRUE)

dev.off()

# Rossi signature

rossi.sig = rowMeans(fingerprint.database.human.15[,as.character(rossi[,1])])
rossi.sig.fingerprint<-rep(0,length(names(rossi.sig)))
names(rossi.sig.fingerprint)<-names(rossi.sig)
rossi.sig.fingerprint[rossi.sig>signature.threshold]<-1
rossi.sig.fingerprint[rossi.sig<(-signature.threshold)]<-(-1)

rossi.sig.difference<-thresholded.fingerprints.full[(abs(rossi.sig)>signature.threshold),] - rossi.sig.fingerprint[(abs(rossi.sig)>signature.threshold)]
rossi.sig.distance<-colSums(abs(rossi.sig.difference))

pdf(file = "Rossi_RiPS_distance_from_Rossi_sig.pdf")
par(mfcol = c(2,1), mar = c(0, 4, 4, 0))
hist(rossi.sig.distance, xlim = c(0,300), breaks = seq(0,300,1), main = "", xlab = "", margin = c(0, 4, 4, 0) )
par(mar = c(7, 4, 4, 0))
hist(rossi.sig.distance[pluripotents], xlim = c(0,300), breaks = seq(0,300,1), main = "", xlab = "distance from RiPS signature", col = "green", margin = c(5, 4, 4, 2))
hist(rossi.sig.distance[as.character(rossi[,1])], xlim = c(0,300), breaks = seq(0,300,1), main = "", xlab = "distance from RiPS signature", col = "red", margin = c(5, 4, 4, 2), add = TRUE)
dev.off()

# extract top arrays 7 of top 10 are the Rossi arrays themselves
as.character(rossi[,1]) %in% names(head(rossi.sig.distance[order(rossi.sig.distance)],10))

#GSM462815_Human induced pluripotent stem cells
#Transcriptional signature and memory retention of human-induced pluripotent stem cells
#PLoS One 2009 Marchetto MC, Yeo GW, Kainohana O, Marsala M, Gage FH, Muotri AR
#
#GSM551195_iPS derived from: PB CD34+ blood
#Loh YH, Hartung O, Li H, Guo C et al. Reprogramming of T cells from human peripheral blood. Cell Stem Cell 2010 Jul 2;7(1):15-9. PMID: 20621044 (G Daley)
#
#GSM553717_Reprogrammed iPSC line hiPSC 3975
#Tchieu J, Plath K
#
#GSM452728_iPS cells derived from GD Nat Genet. 2009
#Doi A, Park I, Wen B, Murakami P, Daley GQ, Feinberg AP
#
#GSM556994_human ES cell line (Hsf1 p51)
#Chin MH, Lowry WE, Plath K

load("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/GEOdescription.RData")

top50<-as.data.frame(unlist(description[names(head(rossi.sig.distance[order(rossi.sig.distance)],50))]))
top50[,2]<-rownames(top50)
colnames(top50)<-c("Description", "GEO reference")
write.table(top50, file = "rossitop50.txt", row.names = FALSE, quote = FALSE, sep = "\t")


rossi.sig.fingerprint.sub<-c(rossi.sig.fingerprint[1:347], rep(FALSE, 491-347)
pdf("rossi.sig.txt")
heatmap(fingerprint.database.human.15[abs(rossi.sig.fingerprint.sub)>0, as.character(rossi[,1])], scale = "none", margin = c(5,25), Colv = NA, cexRow = 0.5, cexCol = 0.5, labCol = rossi[,2])
dev.off()


rossi.pluri<-cbind(fingerprint.database.human.15[(abs(pluripotent.test.signature)>signature.threshold), as.character(rossi[,1])], pluripotent.test.signature.fingerprint[(abs(pluripotent.test.signature)>signature.threshold)])


heatmap(rossi.pluri, scale = "none", margin = c(5,40), Colv = NA, cexRow = 0.5, cexCol = 0.5, labCol = c(as.character(rossi[,2]), "Pluripotent consensus"))
dev.off()

# Figure for paper
# show position of conrad ES cells
conradES<-c("GSM282009", "GSM282010", "GSM282011")
conradGS<-c("GSM282012", "GSM282013", "GSM282008")

pdf("figure for paper")
par(mfcol = c(2,1), mar = c(0, 4, 4, 0))
hist(pluripotent.distance, xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "", margin = c(0, 4, 4, 0) )
par(mar = c(7, 4, 4, 0))
hist(pluripotent.distance.record[pluripotents], xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "distance from pluripotency signature", col = "pink", margin = c(5, 4, 4, 2))
hist(pluripotent.distance.record[oligopotents], xlim = c(0,200), breaks = seq(0,200,1), col = "white", main = "", add = TRUE)
hist(pluripotent.distance.record[conradES], xlim = c(0,200), breaks = seq(0,200,1), col = "green", main = "", add = TRUE)
hist(pluripotent.distance.record[conradGS], xlim = c(0,200), breaks = seq(0,200,1), col = "blue", main = "", add = TRUE)
hist(pluripotent.distance.record[testis], xlim = c(0,200), breaks = seq(0,200,1), col = "red", main = "", add = TRUE)


dev.off()

GEOpossibleES<-names(pluripotent.distance.record)[pluripotent.distance.record<44]

GSE6885<-read.delim("GSE6885.txt", header = FALSE)
GSE6885<-as.character(GSE6885[,1])
GSE6885.BPE<-hist(pluripotent.distance.record[GSE6885[1:5]], xlim = c(0,200), breaks = seq(0,200,1), col = "green", main = "", add = TRUE)
GSE6885.BPLER<-hist(pluripotent.distance.record[GSE6885[6:11]], xlim = c(0,200), breaks = seq(0,200,1), col = "black", main = "", add = TRUE)
GSE6885.HME<-hist(pluripotent.distance.record[GSE6885[c(12, 14:16)]], xlim = c(0,200), breaks = seq(0,200,1), col = "red", main = "", add = TRUE)
GSE6885.HMLER<-hist(pluripotent.distance.record[GSE6885[c(13,17:21)]], xlim = c(0,200), breaks = seq(0,200,1), col = "yellow", main = "", add = TRUE)

# plot histogram on log scale


#par(mfcol = c(2,1), mar = c(0, 4, 4, 0))
temp.full<-hist(pluripotent.distance, xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "", margin = c(0, 4, 4, 0) )
#par(mar = c(7, 4, 4, 0))
temp<-hist(pluripotent.distance.record[pluripotents], xlim = c(0,200), breaks = seq(0,200,1), main = "", xlab = "distance from pluripotency signature", col = "pink", margin = c(5, 4, 4, 2), plot = FALSE)
temp2<-hist(pluripotent.distance.record[conradES], xlim = c(0,200), breaks = seq(0,200,1), col = "green", main = "", add = TRUE, plot = FALSE)
temp3<-hist(pluripotent.distance.record[conradGS], xlim = c(0,200), breaks = seq(0,200,1), col = "blue", main = "", add = TRUE, plot = FALSE)
temp4<-hist(pluripotent.distance.record[testis], xlim = c(0,200), breaks = seq(0,200,1), col = "red", main = "", add = TRUE, plot = FALSE)
pdf("figure for paper_logged.pdf")
par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 

#plot(x = temp$mids, y = 1+log(temp$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "pink", margin = c(5, 4, 4, 2), lwd = 2, ylim = c(0,6), xlim = c(0,200)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)

points(x = temp2$mids, y = 1+log(temp2$counts,2), type = "h", col = "green", lwd = 2.5) 
points(x = temp3$mids, y = 1+log(temp3$counts,2), type = "h", col = "blue", lwd = 2.5) 
points(x = temp4$mids, y = 1+log(temp4$counts,2), type = "h", col = "red", lwd = 2.5) 
dev.off()

# which arrays are within 95% of the maximum values
description[match(names(pluripotent.distance.record)[pluripotent.distance.record<29], names(description))]

nonES<-read.delim("Non stem cell arrays.txt", sep = "-", header = FALSE, skip = 3)

load("/Users/GabrielAltschuler/Documents/Databases/GEOmetadb/gsm_characteristics.RData")

###### analysis of Ovarian Stem cell signature
GSE9891<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Ovarian cancer study/GSE9891.txt", header = FALSE)

as.character(GSE9891[,1]) %in% names(pluripotent.distance.record)
temp.GSE9891<-hist(pluripotent.distance.record[as.character(GSE9891[,1])], breaks = seq(0,200,1), plot = FALSE)


pdf("Ovarian cancer arrays.pdf")
par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)
points(x = temp.GSE9891$mids, y = 1+log(temp.GSE9891$counts,2), type = "h", col = "green", lwd = 2.5)
dev.off()


#### Write list of ES/iPS arrays
write.table(gsm.char[match(pluripotents, gsm.char$gsm),], file = "pluripotent arrays.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

points(x = GSE6885.BPE$mids, y = 1+log(GSE6885.BPE$counts,2), type = "h", col = "green", lwd = 2.5)
points(x = GSE6885.BPLER$mids, y = 1+log(GSE6885.BPLER$counts,2), type = "h", col = "red", lwd = 2.5)
points(x = GSE6885.HME$mids, y = 1+log(GSE6885.HME$counts,2), type = "h", col = "purple", lwd = 2.5)
points(x = GSE6885.HMLER$mids, y = 1+log(GSE6885.HMLER$counts,2), type = "h", col = "yellow", lwd = 2.5)

pluripotents.description<-gsm.char[match(pluripotents,gsm.char[,1]),]

write.table(pluripotents.description, file = "pluripotent.table.txt", sep = "\t", quote = FALSE, row.names = FALSE)

###### analysis of Cancer Stem cell Arrays
GOarrays<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/GO presentation/GO grant records.txt", stringsAsFactors = FALSE)
GOarrays.valid<-GOarrays[(GOarrays$GSM %in% names(pluripotent.distance.record)),]
GOarrays.CancerStem<-GOarrays.valid[grep("CancerStem", GOarrays.valid$Status),]
GOarrays.Cancer<-GOarrays.valid[grep("^Cancer$", GOarrays.valid$Status),]


temp.GOarrays.CancerStem<-hist(pluripotent.distance.record[as.character(GOarrays.CancerStem$GSM)], breaks = seq(0,200,1), plot = FALSE)
temp.GOarrays.Cancer<-hist(pluripotent.distance.record[as.character(GOarrays.Cancer$GSM)], breaks = seq(0,200,1), plot = FALSE)


# pdf("cancer stem cell arrays.pdf")
par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)
points(x = temp.GOarrays.CancerStem$mids, y = 1+log(temp.GOarrays.CancerStem$counts,2), type = "h", col = "green", lwd = 2.5)
points(x = temp.GOarrays.Cancer$mids, y = 1+log(temp.GOarrays.Cancer$counts,2), type = "h", col = "red", lwd = 2.5)
dev.off()

######
# preparing stepped through figure for PQG talk
pdf("figure for PQG.pdf")
par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 

par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)

par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)
points(x = temp2$mids, y = 1+log(temp2$counts,2), type = "h", col = "green", lwd = 2.5) 

par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)
points(x = temp2$mids, y = 1+log(temp2$counts,2), type = "h", col = "green", lwd = 2.5) 
points(x = temp3$mids, y = 1+log(temp3$counts,2), type = "h", col = "blue", lwd = 2.5) 

par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)
points(x = temp2$mids, y = 1+log(temp2$counts,2), type = "h", col = "green", lwd = 2.5) 
points(x = temp3$mids, y = 1+log(temp3$counts,2), type = "h", col = "blue", lwd = 2.5) 
points(x = temp4$mids, y = 1+log(temp4$counts,2), type = "h", col = "red", lwd = 2.5) 
dev.off()

######
# preparing stepped through figure for CD133 for Rama
cd133pos.hist<-hist(pluripotent.distance.record[cd133.pos], breaks = seq(0,200,1), plot = FALSE)
cd133neg.hist<-hist(pluripotent.distance.record[cd133.neg], breaks = seq(0,200,1), plot = FALSE)

pdf("CD133 vs pluripotency figure.pdf")
par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 

par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)

par(mar = c(20, 4, 4, 4))
plot(x = temp.full$mids, y = 1+log(temp.full$counts,2), type = "h", main = "", xlab = "distance from pluripotency signature", col = "grey", lwd = 2.5, xlim = c(0,200), ylim = c(0,12.5), ylab = "1+log(frequency)", yaxp = c(1,12,11)) 
points(x = temp$mids, y = 1+log(temp$counts,2), type = "h", col = "pink", lwd = 2.5)
points(x = cd133pos.hist$mids, y = 1+log(cd133pos.hist$counts,2), type = "h", col = "green", lwd = 2.5) 
points(x = cd133neg.hist$mids, y = 1+log(cd133neg.hist$counts,2), type = "h", col = "blue", lwd = 2.5) 
dev.off()

save.image("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Data/brief communication dataset.RData")

