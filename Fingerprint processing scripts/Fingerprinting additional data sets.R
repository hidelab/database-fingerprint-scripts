# preparing genelists for fingerprining run of extra gene lists
source("/Users/GabrielAltschuler/Documents/Databases/gabriel functions.R")
library(annotationTools)
homologeneFile <- '/Users/GabrielAltschuler/Documents/Databases/Homologene/homologene.data'
homologene <- read.delim(homologeneFile, header=F)


# GeneSigDb delayed for now
# 2x Aedin's lists
# Melanoma lists
# Fabio's RNAi hit list
# Peggy's list
# miR34 intersection of PD/DOWN for both cell lines
# Langmoen GBM up/down lists


OvarianDifferentiation<-as.character(read.delim("/Users/GabrielAltschuler/Documents/Projects/Ovarian cancer study/OvarianDifferentiationEntrez.txt", header = FALSE)[,1])
OvarianStemLike<-as.character(read.delim("/Users/GabrielAltschuler/Documents/Projects/Ovarian cancer study/OvarianStemLikeEntrez.txt", header = FALSE)[,1])

CD133_common_Neg_peaks.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Melanoma/Gene lists/Common_Neg_Differential_peaks.txt")
CD133_common_Pos_peaks.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Melanoma/Gene lists/Common_Pos_Differential_peaks.txt")
CD133_common_Neg_peaks<-genename2Entrez(CD133_common_Neg_peaks.frame$ID, species = "human")
CD133_common_Pos_peaks<-genename2Entrez(CD133_common_Pos_peaks.frame$ID, species = "human")
CD133_common_Neg_peaks<-CD133_common_Neg_peaks[!(is.na(CD133_common_Neg_peaks))]
CD133_common_Pos_peaks<-CD133_common_Pos_peaks[!(is.na(CD133_common_Pos_peaks))]

siRNABreastCancerValidatedHits.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Lieberman/RNAi screen/validated_154_siRNAs.txt")
siRNABreastCancerHits.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Lieberman/RNAi screen/Petrocca_siRNA_hits.txt")
siRNABreastCancerValidatedHits<-siRNABreastCancerValidatedHits.frame$Entrez.Gene.ID
siRNABreastCancerHits<-siRNABreastCancerHits.frame$Entrez.Gene.ID

CottonWorkersMouseModel.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/CottonWorkers/GeneSig101_entrez.txt")
CottonWorkersMouseModel<-CottonWorkersMouseModel.frame$EntrezGene.ID
CottonWorkersMouseModel<-entrezUpdate(CottonWorkersMouseModel)
CottonWorkersMouseModel.Hs<-getHOMOLOG(CottonWorkersMouseModel, 9606, homologene)
# take first value
for (i in 1:length(CottonWorkersMouseModel.Hs)){
	CottonWorkersMouseModel.Hs[i]<-unlist(CottonWorkersMouseModel.Hs[i])[1]		}
CottonWorkersMouseModel.Hs<-unlist(CottonWorkersMouseModel.Hs)
CottonWorkersMouseModel.Hs<-CottonWorkersMouseModel.Hs[!(is.na(CottonWorkersMouseModel.Hs))]


PD_HCT116_AND_K562.final.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Lieberman/miR34/Files for final presentations/Data/HCT116 and K562 PD_2.5.txt", sep = " ")
DOWN_HCT116_AND_K562.final.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Lieberman/miR34/Files for final presentations/Data/HCT116 and K562 down.txt", sep = " ")
miR34_PD_DOWN_Intersection<-intersect(PD_HCT116_AND_K562.final.frame$EntrezID, DOWN_HCT116_AND_K562.final.frame$EntrezID)

LangmoenGBM_Up.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Langmoen/Data/results/upregulated_1PFP.txt")
LangmoenGBM_Down.frame<-read.delim("/Users/GabrielAltschuler/Documents/Projects/Langmoen/Data/results/downregulated_1PFP.txt")
LangmoenGBM_Up<-LangmoenGBM_Up.frame$EntrezID
LangmoenGBM_Down<-LangmoenGBM_Down.frame$EntrezID

ExperimentalGeneSigs.Hs.gs<-list(
							OvarianDifferentiation = entrezUpdate(OvarianDifferentiation),
							OvarianStemLike = entrezUpdate(OvarianStemLike),
							CD133_common_Neg_peaks = entrezUpdate(CD133_common_Neg_peaks),
							CD133_common_Pos_peaks = entrezUpdate(CD133_common_Pos_peaks),
							siRNABreastCancerValidatedHits = entrezUpdate(siRNABreastCancerValidatedHits),
							siRNABreastCancerHits = entrezUpdate(siRNABreastCancerHits),
							CottonWorkersMouseModel = entrezUpdate(CottonWorkersMouseModel.Hs),
							miR34_PD_DOWN_Intersection = entrezUpdate(miR34_PD_DOWN_Intersection),
							LangmoenGBM_Up = entrezUpdate(LangmoenGBM_Up),
							LangmoenGBM_Down = entrezUpdate(LangmoenGBM_Down)
							)

# Converting for mouse							
allHsEntrezID<-as.data.frame(unique(unlist(ExperimentalGeneSigs.Hs.gs)))
allHsEntrezID[,1]<-as.character(allHsEntrezID[,1])
humanGenes <- allHsEntrezID[,1]
# Taxonomy info from http://www.ncbi.nlm.nih.gov/Taxonomy/
mouseHomologs <- getHOMOLOG(humanGenes, 10090, homologene)
# take just the first id (tends not to be the riken or LOC number)
mouseHomologs.unique <- mouseHomologs

for (i in 1:length(mouseHomologs.unique)){
	mouseHomologs.unique[i]<-unlist(mouseHomologs.unique[i])[1]		}

allHsEntrezID[,2]<-unlist(mouseHomologs.unique)
colnames(allHsEntrezID)<-c("human", "mouse")

# create conversion frame
allHsEntrezID.frame<-allHsEntrezID
rownames(allHsEntrezID.frame)<-allHsEntrezID.frame[,1]
allHsEntrezID.frame[,2]<-entrezUpdate(allHsEntrezID.frame[,2])

# proportion NA is 6.6%
(sum(is.na(allHsEntrezID[,2]))/length(allHsEntrezID[,2]))*100

# now use this table to convert gene ids

ExperimentalGeneSigs.Mm.gs<-vector("list", 0)
for (i in 1:length(ExperimentalGeneSigs.Hs.gs)){	
	temp<-allHsEntrezID.frame[unlist(ExperimentalGeneSigs.Hs.gs[i]),2]
	ExperimentalGeneSigs.Mm.gs <- append(ExperimentalGeneSigs.Mm.gs, list(unique(na.omit(unlist(as.character(temp))))))
		}

names(ExperimentalGeneSigs.Mm.gs)<-names(ExperimentalGeneSigs.Hs.gs)

# put Peggy's list back in the original form

ExperimentalGeneSigs.Mm.gs$CottonWorkersMouseModel<-entrezUpdate(CottonWorkersMouseModel)

# save files

save(ExperimentalGeneSigs.Hs.gs, file = "/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Hs.gs")
save(ExperimentalGeneSigs.Mm.gs, file = "/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Mm.gs")

###############################
# updated 12th April 2010
# Adding 3 lists from Rama and the Langmoen larger up/down signature
source("/Users/GabrielAltschuler/Dropbox/fingerprint/scripts/gabriel functions.R")
load("/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Hs.gs")
load("/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Mm.gs")

# save old sets as backup
save(ExperimentalGeneSigs.Hs.gs, file = "/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Hs.gs.old")
save(ExperimentalGeneSigs.Mm.gs, file = "/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Mm.gs.old")

# load new genes
# Langmoen
glioma.G<-scan("/Users/GabrielAltschuler/Documents/Projects/Langmoen/Survival Analysis/data/GeneSets/Glioma_G.txt", "character", skip = 1)
glioma.H<-scan("/Users/GabrielAltschuler/Documents/Projects/Langmoen/Survival Analysis/data/GeneSets/Glioma_H.txt", "character", skip = 1)
glioma.G.entrez<-genename2Entrez(glioma.G, "human")
glioma.H.entrez<-genename2Entrez(glioma.H, "human")
glioma.G.entrez<-as.character(na.omit(glioma.G.entrez))
glioma.H.entrez<-as.character(na.omit(glioma.H.entrez))

# Melanoma lists
cd133hi_downreg<-readLines("/Users/GabrielAltschuler/Documents/Projects/Melanoma/Data/genelists/cd133hi_downreg.txt")
cd133hi_upreg<-readLines("/Users/GabrielAltschuler/Documents/Projects/Melanoma/Data/genelists/cd133hi_upreg.txt")


cd133hi_downreg.entrez<-genename2Entrez(cd133hi_downreg, "human")
cd133hi_upreg.entrez<-genename2Entrez(cd133hi_upreg, "human")

# add missing ID
cd133hi_downreg.entrez<-c(cd133hi_downreg.entrez, 57730)

nyeso1_ip_shortlist_ez<-readLines("/Users/GabrielAltschuler/Documents/Projects/Melanoma/Data/genelists/nyeso1_ip_shorlist_ez.txt")

newSets.Hs<-list(
				glioma.G = entrezUpdate(glioma.G.entrez),
				glioma.H = entrezUpdate(glioma.H.entrez),
				cd133hi_downreg = entrezUpdate(cd133hi_downreg.entrez),
				cd133hi_upreg = entrezUpdate(cd133hi_upreg.entrez),
				nyeso1_ip_shortlist_ez = entrezUpdate(nyeso1_ip_shortlist_ez)
							)

# remove NAs
charNAremove<-function(x){x<-x[!(x == "NA")]}
newSets.Hs<-lapply(newSets.Hs, charNAremove)
newSets.Hs<-lapply(newSets.Hs, na.omit)
newSets.Hs<-lapply(newSets.Hs, as.character)

library(annotationTools)
homologeneFile <- '/Users/GabrielAltschuler/Documents/Databases/Homologene/homologene.data'
homologene <- read.delim(homologeneFile, header=F)




# Converting for mouse							
allHsEntrezID<-as.data.frame(unique(unlist(newSets.Hs)))
allHsEntrezID[,1]<-as.character(allHsEntrezID[,1])
humanGenes <- allHsEntrezID[,1]
# Taxonomy info from http://www.ncbi.nlm.nih.gov/Taxonomy/
mouseHomologs <- getHOMOLOG(humanGenes, 10090, homologene)
# take just the first id (tends not to be the riken or LOC number)
mouseHomologs.unique <- mouseHomologs

for (i in 1:length(mouseHomologs.unique)){
	mouseHomologs.unique[i]<-unlist(mouseHomologs.unique[i])[1]		}

allHsEntrezID[,2]<-unlist(mouseHomologs.unique)
colnames(allHsEntrezID)<-c("human", "mouse")

# create conversion frame
allHsEntrezID.frame<-allHsEntrezID
rownames(allHsEntrezID.frame)<-allHsEntrezID.frame[,1]
allHsEntrezID.frame[,2]<-entrezUpdate(allHsEntrezID.frame[,2])

# proportion NA is 12.9%
(sum(is.na(allHsEntrezID[,2]))/length(allHsEntrezID[,2]))*100

# now use this table to convert gene ids

newSets.Mm<-vector("list", 0)
for (i in 1:length(newSets.Hs)){	
	temp<-allHsEntrezID.frame[unlist(newSets.Hs[i]),2]
	newSets.Mm <- append(newSets.Mm, list(unique(entrezUpdate(na.omit(unlist(as.character(temp)))))))
		}

names(newSets.Mm)<-names(newSets.Hs)

ExperimentalGeneSigs.Hs.gs<-append(ExperimentalGeneSigs.Hs.gs[1:10], newSets.Hs)
ExperimentalGeneSigs.Mm.gs<-append(ExperimentalGeneSigs.Mm.gs[1:10], newSets.Mm)


# bug fix - need to remove "NA" as well as NA
charNAremove<-function(x){x<-x[!(x == "NA")]}
ExperimentalGeneSigs.Hs.gs <- lapply(ExperimentalGeneSigs.Hs.gs, charNAremove)
ExperimentalGeneSigs.Mm.gs <- lapply(ExperimentalGeneSigs.Mm.gs, charNAremove)



# save locally
save(ExperimentalGeneSigs.Hs.gs, file = "/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Hs.gs")
save(ExperimentalGeneSigs.Mm.gs, file = "/Users/GabrielAltschuler/Documents/Databases/Gene sets/ExperimentalGeneSigs.Mm.gs")

# save to repository
save(ExperimentalGeneSigs.Hs.gs, file = "/Users/GabrielAltschuler/Dropbox/fingerprint/data/ExperimentalGeneSigs.Hs.gs")
save(ExperimentalGeneSigs.Mm.gs, file = "/Users/GabrielAltschuler/Dropbox/fingerprint/data/ExperimentalGeneSigs.Mm.gs")