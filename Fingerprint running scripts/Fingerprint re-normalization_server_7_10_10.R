# re-running normalization on server

# loading and analysis of fingerprints
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/")

fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints"
setwd(fingerpath)
source("/home/galtschu2/Documents/Databases/gabriel functions.R")

if (!(exists("kegg_wiki_TR_static.Hs.gs"))){
	load("/home/galtschu2/Documents/Databases/Gene sets/kegg_wiki_TR_static.Hs.gs", .GlobalEnv)
	}
if (!(exists("kegg_wiki_TR_static.Mm.gs"))){
	load("/home/galtschu2/Documents/Databases/Gene sets/kegg_wiki_TR_static.Mm.gs", .GlobalEnv)
	}

if (!(exists("chipframe"))){
	load("/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe.R", .GlobalEnv)
	}


# only load complete fingerprints (111348) this takes some time (10mins) on the server
# with update now have 115813 files

files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
files<-files[grep("KEGG_Wiki_static_", files)]

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
	
names(SCG)<-gsub("KEGG_Wiki_static_", "", files)
SCG.frame<-t(as.data.frame(SCG))
colnames(SCG.frame)<-names(kegg_wiki_TR_static.Hs.gs)


# re-analysis of SCG

load("SCG.frame.RData")

# now need to normalize this according to the chip data


platform<-vector("list", length(frame))
pb <- txtProgressBar(min = 0, max = length((frame)), style = 3)
for (i in 1:length(frame)){
	platform[i]<-frame[[i]]$platform
	setTxtProgressBar(pb, i)
	}

names(platform)<-gsub("KEGG_Wiki_static_", "", files)
platform.frame<-data.frame(GEO = names(platform), Platform = unlist(Platform))
save(platform.frame, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/platform.frame.RData")


chiplengths <- c(16477, 12170, 5779, 18110, 12056, 7738, 8341, 12170, 12056, 7030) 
names(chiplengths)<-c("GPL1261", "GPL339", "GPL340", "GPL570", "GPL571", "GPL81", "GPL8300", "GPL8321", "GPL96", "GPL97")
species <- c("mouse", "mouse", "mouse", "human", "human", "mouse", "human", "mouse", "human", "human")
names(species)<-names(chiplengths)



platform.lengths <- vector("list", length(platform))
pb <- txtProgressBar(min = 0, max = length((frame)), style = 3)
for (i in 1:length(platform.lengths)){
	platform.lengths[i]<-chiplengths[platform[[i]]]
	setTxtProgressBar(pb, i)
	}

# need to re-calculate the overlaps

chipgenes<-vector("list", length(chiplengths))
names(chipgenes)<-names(chiplengths)
for (i in 1:length(chipgenes)){
	chipgenes[i]<-list(unique(gsub("_at", "", chipframe[[names(chiplengths)[i]]]$ann$EntrezID)))
	}


kegg_wiki.Hs.overlaplength<-function(x){
	temp<-vector("numeric", length(kegg_wiki_TR_static.Hs.gs))
	for (i in 1:length(kegg_wiki_TR_static.Hs.gs)){
		temp[i]<-length(intersect(x, kegg_wiki_TR_static.Hs.gs[[i]]))		}
	return(temp)
	}

kegg_wiki.Mm.overlaplength<-function(x){
	temp<-vector("numeric", length(kegg_wiki_TR_static.Mm.gs))
	for (i in 1:length(kegg_wiki_TR_static.Mm.gs)){
		temp[i]<-length(intersect(x, kegg_wiki_TR_static.Mm.gs[[i]]))		}
	return(temp)
	}


overlaps<-vector("list", length(chiplengths))
names(overlaps)<-names(chiplengths)
for (i in 1:length(chipgenes)){
	if (species[names(overlaps[i])] == "human"){
		overlaps[i]<-list(kegg_wiki.Hs.overlaplength(chipgenes[[names(chiplengths)[i]]]))
		}
	if (species[names(overlaps[i])] == "mouse"){
		overlaps[i]<-list(kegg_wiki.Mm.overlaplength(chipgenes[[names(chiplengths)[i]]]))
		} 
	}


# Problem - shouldn't really use pathways for which the overlap is <5 on all chips...

table(overlaps[[1]])

names(kegg_wiki_TR_static.Mm.gs)[overlaps[[3]] %in% 0]

# it may be that the best thing to do is to normalize each chip type separately - then you have the same representation of genes with each chip. However, can you still compare properly across platform...? Hmmm interesting...
# need to set scroe to NA if overlap is <5 but this 

# are there some pathways for which there is never any significant overlap...?
overlaps.frame<-as.data.frame(overlaps)

names(kegg_wiki_TR_static.Mm.gs)[rowSums(overlaps.frame<5) == 10]
# [1] "00300 Lysine biosynthesis"                   
# [2] "00401 Novobiocin biosynthesis"               
# [3] "00471 D-Glutamine and D-glutamate metabolism"
# [4] "00472 D-Arginine and D-ornithine metabolism" 
# [5] "00643 Styrene degradation"                   
# [6] "00780 Biotin metabolism"                     
# [7] "00785 Lipoic acid metabolism"                
# [8] "amino acid conjugation of benzoic acid"      
# [9] "Polyol pathway"                              
#[10] "Serotonin Receptor 2 -> STAT3 signaling"     
#[11] "Alpha6 Beta4 Integrin DOWN"                  
#[12] "Hedgehog DOWN"                               
#[13] "Kit Receptor DOWN"                           
#[14] "IL-9 DOWN"                                   

#These are mainly pathways with only 4 members originally, apary from Alpha 6Beta4 Integrin DOWN which perhaps needs an entrezUpdate - actually this is not the case but it includes predicted genes some of which may not be on the affy chips.
# There are more pathways that are not represented by at least 5 members in at least 8 of the chip types
names(kegg_wiki_TR_static.Mm.gs)[rowSums(overlaps.frame<5) > 2]
# [1] "00061 Fatty acid biosynthesis"                            
# [2] "00072 Synthesis and degradation of ketone bodies"         
# [3] "00300 Lysine biosynthesis"                                
# [4] "00400 Phenylalanine, tyrosine and tryptophan biosynthesis"
# [5] "00401 Novobiocin biosynthesis"                            
# [6] "00430 Taurine and hypotaurine metabolism"                 
# [7] "00440 Aminophosphonate metabolism"                        
# [8] "00460 Cyanoamino acid metabolism"                         
# [9] "00471 D-Glutamine and D-glutamate metabolism"             
#[10] "00472 D-Arginine and D-ornithine metabolism"              
#[11] "00625 Tetrachloroethene degradation"                      
#[12] "00643 Styrene degradation"                                
#[13] "00730 Thiamine metabolism"                                
#[14] "00750 Vitamin B6 metabolism"                              
#[15] "00780 Biotin metabolism"                                  
#[16] "00785 Lipoic acid metabolism"                             
#[17] "00791 Atrazine degradation"                               
#[18] "00830 Retinol metabolism"                                 
#[19] "00930 Caprolactam degradation"                            
#[20] "00940 Phenylpropanoid biosynthesis"                       
#[21] "Synthesis and Degradation of Ketone Bodies"               
#[22] "amino acid conjugation of benzoic acid"                   
#[23] "ACE Inhibitor Pathway"                                    
#[24] "Arachidonate Epoxygenase / Epoxide Hydrolase"             
#[25] "Catalytic cycle of mammalian FMOs"                        
#[26] "Polyol pathway"                                           
#[27] "Arylamine metabolism"                                     
#[28] "Serotonin Receptor 2 -> STAT3 signaling"                  
#[29] "Alpha6 Beta4 Integrin DOWN"                               
#[30] "Hedgehog DOWN"                                            
#[31] "Kit Receptor DOWN"                                        
#[32] "IL-9 DOWN"                                                
#[33] "{ACY1,11}"                                                
#[34] "{FLI1,10}"                                                
#[35] "{GNAT3,28}"                                               
#[36] "{HSP90B1,11}"                                             
#[37] "{KCNB1,10}"                                               
#[38] "{PAK2,10}"                                                
#[39] "{RAD21,11}"    

# next step is to use pnorm to convert all of the ESs to NESs based on the pnorm distribution for that chip and pathway overlap. With an additional caveat of NA score if pathway has <5 members in the overlap.

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

# These parameters were actually saved using a previous session
fingerprintParameters<-list(species = species, overlaps = overlaps, chiplengths = chiplengths, par.mean = par.mean, par.sd.b = par.sd.b, par.sd.power = par.sd.power, par.sd.a = par.sd.a, platform7_10_10 = platform)

save(fingerprintParameters, file = "/home/galtschu2/Documents/Projects/Fingerprinting/fingerprintParameters.RData")

# also saving overlaps.frame separately
save(overlaps.frame, file = "/home/galtschu2/Documents/Projects/Fingerprinting/overlaps.frame.RData")

# and for a genelist with intersect OVERLAP, normalized score for SCORE is
#pnorm(SCORE, mean = par.mean*(1-(OVERLAP^-1)), sd = par.sd.a+(par.sd.b*(OVERLAP^par.sd.power))))


Norm.SCG<-vector("list", length(SCG))
pb <- txtProgressBar(min = 0, max = length(SCG), style = 3)
for (i in 1:length(SCG)){
	intersect <- overlaps[[platform[[i]]]]
	intersect[intersect<5] <- NA
	mean<-par.mean[platform[[i]]]*(1-(intersect^-1))
	sd<- par.sd.a+(par.sd.b[platform[[i]]]*(intersect^par.sd.power[platform[[i]]]))
	Norm.SCG[i]<-list(pnorm(SCG[[i]], mean = mean, sd = sd))
	setTxtProgressBar(pb, i)
	}

names(Norm.SCG)<-gsub("KEGG_Wiki_static_", "", files)

Norm.SCG.frame<-t(as.data.frame(Norm.SCG))
colnames(Norm.SCG.frame)<-names(kegg_wiki_TR_static.Hs.gs)

# now the normalized matrix had been obtained need to re-normalize this using the ecdf to obtain the norm-norm matrix.
Norm.SCG.frame<-t(Norm.SCG.frame)

# need to produce the ECDF, split into 14 parts to run parallel as takes 15 seconds each
Norm.SCG.matrix<-as.matrix((Norm.SCG.frame))

#problem is that pathways with all NAs cannot be processed into an ecdf, therefore need to remove or replace with dummy numbers

Norm.SCG.matrix[names(kegg_wiki_TR_static.Mm.gs)[rowSums(overlaps.frame<5) == 10],]<-1

Norm.SCG.ecdf<-apply(Norm.SCG.matrix, 1, ecdf)


save(Norm.SCG.matrix, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/Norm.SCG.matrix7_10_10.Rdata")
save(Norm.SCG.ecdf, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/Norm.SCG.ecdf7_10_10.Rdata")

Norm.SCG.Norm.matrix<-matrix(ncol = length(Norm.SCG.matrix[1,]), nrow = length(Norm.SCG.matrix[,1]))
colnames(Norm.SCG.Norm.matrix)<-colnames(Norm.SCG.matrix)
rownames(Norm.SCG.Norm.matrix)<-rownames(Norm.SCG.matrix)

for (i in 1:length(Norm.SCG.ecdf)){
	Norm.SCG.Norm.matrix[i,]<-Norm.SCG.ecdf[[i]](Norm.SCG.matrix[i,])
	}

Norm.SCG.Norm.matrix[names(kegg_wiki_TR_static.Mm.gs)[rowSums(overlaps.frame<5) == 10],]<-NA
	
save(Norm.SCG.Norm.matrix, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/Norm.SCG.Norm.matrix7_10_10.Rdata")

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data")
pdf(file = "Norm_SCG_Norm_hist7_10_10.pdf")
par(mfcol = c(2,1))

breaks<-0.005*(0:200)
for (i in 1:length(Norm.SCG.matrix[,1])){
	hist(Norm.SCG.matrix[i,], xlim = c(0,1), breaks = breaks, main = paste(rownames(Norm.SCG.matrix)[i], "Norm_SCG", sep = " "))
	hist(Norm.SCG.Norm.matrix[i,], xlim = c(0,1), breaks = breaks, main = paste(rownames(Norm.SCG.Norm.matrix)[i], "Norm_SCG_Norm", sep = " "))
	}
dev.off()	

# these plots look identical to those produced during the last normalization so assume okay
# the re-normalization appears to have worked in terms of creating a uniform distribution of enrichment scores for each pathway. The challenge now is to cluster effectively.

zeropathways<-names(kegg_wiki_TR_static.Mm.gs)[rowSums(overlaps.frame<5) == 10]


# free up memory

memory<-vector("numeric", length(objects))
objects <- ls()
for (i in 1:length(objects)){
	memory[i] <- object.size(get(objects[i]))
	}

names(memory)<-objects

head(memory[order(memory, decreasing = TRUE)])


save(Norm.SCG.frame, file = "Norm.SCG.frame7_10_10.RData")

rm(Norm.SCG, Norm.SCG.frame, frame, SCG.frame, SCG)

# now need to re-normalize by chip

# Re-normalizing fingerprints separately for each chip
# problem is that for normalizing each chip internally there is not enough pathway representation in the smaller chips to give any proper values

colSums(overlaps.frame<5)
# GPL340 and GPL97 have over 200 pathways with representation of <5 members, GPL81 is the next highest with 37, these are predominantly metabolic pathways. GPL81 is MGU74Av2, a chip represented in stembase so better to include it in the fingerprint.
# However, exclude GPL340 and GPL97 for now.
# These pathways should be excluded

names(kegg_wiki_TR_static.Mm.gs)[rowSums(overlaps.frame[,c(1,2,4,5,6,7,8,9)]<5) > 0]

#set all of these pathways to 1

Norm.SCG.matrix[names(kegg_wiki_TR_static.Mm.gs)[rowSums(overlaps.frame[,c(1,2,4,5,6,7,8,9)]<5) > 0],]<-1

# now separate out each chip and produce the ecdf - N.B. the save step of the ecdfs was carried out in a different session (cd133)


for (i in 1:length(chiplengths)){
	temp<-Norm.SCG.matrix[,grep(names(chiplengths)[i], platform)]
	assign(paste(names(chiplengths)[i], ".Norm.SCG.matrix", sep = ""), temp)
	}
	
for (i in c(1,2,4,5,6,7,8,9)){
	temp<-apply(get(paste(names(chiplengths)[i], ".Norm.SCG.matrix", sep = "")), 1, ecdf)
	assign(paste(names(chiplengths)[i], ".Norm.SCG.ecdf", sep = ""), temp)
	save(temp, file = paste("/home/galtschu2/Documents/Projects/Fingerprinting/data/", names(chiplengths)[i], ".7_10_10_Norm.SCG.ecdf.RData", sep = ""))
	}

pdf(file = "Norm_hist_by_chip_7_10_10.pdf")	
par(mfcol = c(4,2))
for (j in 1:length(Norm.SCG.matrix[,1])){
	for (i in c(1,2,4,5,6,7,8,9)){	
		hist(get(paste(names(chiplengths)[i], ".Norm.SCG.matrix", sep = ""))[j,], xlim = c(0,1), breaks = breaks, main = paste(names(chiplengths)[i], rownames(Norm.SCG.matrix)[j], sep = "_"), xlab = "Q1 normalized score")
		}
	}
	
dev.off()
# looks the same as previous normalization

# now normalize each chip separately
for (i in c(1,2,4,5,6,7,8,9)){
	temp1<-get(paste(names(chiplengths)[i], ".Norm.SCG.matrix", sep = ""))
	temp<- matrix(ncol = length(temp1[1,]), nrow = length(temp1[,1]))
	colnames(temp)<-colnames(temp1)
	rownames(temp)<-rownames(temp1)
	temp2<-get(paste(names(chiplengths)[i], ".Norm.SCG.ecdf", sep = ""))
	for (j in 1:length(temp2)){
		temp[j,]<-temp2[[j]](temp1[j,])
		}			
	temp[names(kegg_wiki_TR_static.Mm.gs)[rowSums(overlaps.frame[,c(1,2,4,5,6,7,8,9)]<5) > 0],]<-NA
	assign(paste(names(chiplengths)[i], ".N.Norm.SCG", sep = ""), temp)
	}

# re-plotting should all be uniform distributions (except the NAs) - it worked
pdf(file = "Norm_Norm_hist_by_chip7_10_10.pdf")	
par(mfcol = c(4,2))
for (j in 1:length(Norm.SCG.matrix[,1])){
	for (i in c(1,2,4,5,6,7,8,9)){	
		hist(get(paste(names(chiplengths)[i], ".N.Norm.SCG", sep = ""))[j,], xlim = c(0,1), breaks = breaks, main = paste(names(chiplengths)[i], rownames(Norm.SCG.matrix)[j], sep = "_"), xlab = "Q1 normalized score")
		}
	}
	
dev.off()

# re-assemble matrix

ReNorm.matrix<-cbind(
					get(paste(names(chiplengths)[1], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[2], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[4], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[5], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[6], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[7], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[8], ".N.Norm.SCG", sep = "")),
					get(paste(names(chiplengths)[9], ".N.Norm.SCG", sep = ""))					)


# comparison to last time shows only slight changes, as expected from adding about 5000 arrays

save(ReNorm.matrix, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix7_10_10.RData")

##----##
# now quit and restart to free up memory on the server

load("/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix7_10_10.RData")

# re-apply thresholds, higher cutoff as within smaller groups
high<-0.95
low<-0.05
ReNorm.matrix.ternary<-((ReNorm.matrix)>high)-((ReNorm.matrix)<low)
ReNorm.matrix.ternary[is.na(ReNorm.matrix.ternary)]<-0

# create lower threshold dataset ############### run this part in NY ###################
# re-apply thresholds, higher cutoff as within smaller groups
high.10<-0.9
low.10<-0.1
ReNorm.matrix.ternary.10<-((ReNorm.matrix)>high.10)-((ReNorm.matrix)<low.10)
ReNorm.matrix.ternary.10[is.na(ReNorm.matrix.ternary.10)]<-0

high.15<-0.85
low.15<-0.15
ReNorm.matrix.ternary.15<-((ReNorm.matrix)>high.15)-((ReNorm.matrix)<low.15)
ReNorm.matrix.ternary.15[is.na(ReNorm.matrix.ternary.15)]<-0

high.20<-0.80
low.20<-0.20
ReNorm.matrix.ternary.20<-((ReNorm.matrix)>high.20)-((ReNorm.matrix)<low.20)
ReNorm.matrix.ternary.20[is.na(ReNorm.matrix.ternary.20)]<-0

high.25<-0.75
low.25<-0.25
ReNorm.matrix.ternary.25<-((ReNorm.matrix)>high.25)-((ReNorm.matrix)<low.25)
ReNorm.matrix.ternary.25[is.na(ReNorm.matrix.ternary.25)]<-0


########----########
# save ternary matrix, this is an 8.3MB file so can be transferred and sftp'd over to be locally accessed
save(ReNorm.matrix.ternary, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix.ternary7_10_10.RData")
save(ReNorm.matrix.ternary.10, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix.10ternary_7_10_10.RData")
save(ReNorm.matrix.ternary.15, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix.15ternary_7_10_10.RData")
save(ReNorm.matrix.ternary.20, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix.20ternary_7_10_10.RData")
save(ReNorm.matrix.ternary.25, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/ReNorm.matrix.25ternary_7_10_10.RData")


###---###
# Oliver's suggestions 1) Manhatthan clustering 2) alternative distribution to fit for second normalization step. 3) Perhaps try normalizing within each chip set. 4) Use correlation structure within the pathways to construct a Mahonoblis distance matrix and cluster with this as a parameter (how...?)

# need to re-normalize within each chip type rather than across all chips

# see script "re-normalization of fingerprint by chip"

# try a PCA of the entire set - might get a species effect dominating...?



