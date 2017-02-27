# Script to take the 'real' genesets and 
# find the SCG scores for N random gene permutations
# for each chip

# For future reference this should be distributed properly over the CPUs
# Although, would involve a large amount of mem usage so might not be a good idea


source("/home/galtschu2/Documents/Databases/gabriel functions.R")
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data")

# load chip data
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update.RData")

# enter chip types

species <- c("mouse", "mouse", "mouse", "human", "human", "mouse", "human", "mouse", "human", "human", "human")
names(species)<-c("GPL1261", "GPL339", "GPL340", "GPL570", "GPL571", "GPL81", "GPL8300", "GPL8321", "GPL96", "GPL97", "GPL2986")

# load platform types - perhaps better to update this to include the reactome set
load(file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/platform.frame.RData")
platforms<-as.character(unique(platform.frame$Platform))

# load genesets
load("/home/galtschu2/Documents/Databases/Gene sets/kegg_wiki_TR_static.Hs.gs")
load("/home/galtschu2/Documents/Databases/Gene sets/kegg_wiki_TR_static.Mm.gs")


# For next round use reactome.kegg.wiki.netpath human and mouse
# load("/home/galtschu2/Documents/Databases/Gene sets/reactome.kegg.wiki.netpath.static.gs")

permutations <- 10000
# permutations <- 3 # initial test

# create dataframe of random exprs
# perhaps quicker to use SCG to process complete frame
# rather than processing each column indiviually

# also re-started, saving each matrix separately to prevent too much mem usage
# background.SCG <- vector("list", 0)
for (i in 1:length(chipframe))
  {
  # extract list of genes
    chip<-names(chipframe[i])
    print(paste("Processing", chip))
    chipdata<-chipframe[[i]]$ann
    chipgenes<-unique(sub("_at", "", chipdata$EntrezID)) 
  # produce random array, dim (Ngenes x Npermutations)
    randomexprs<-apply(matrix(rep(1:length(chipgenes),permutations),
                              ncol = permutations
                              ),
                      2, sample
                      )
    rownames(randomexprs)<-chipgenes
  # select geneset depending on species
  if (species[chip] == "human")
       {geneset <- kegg_wiki_TR_static.Hs.gs}
  if (species[chip] == "mouse")
        {geneset <- kegg_wiki_TR_static.Mm.gs}

  # run SCG routine  
  	SCG<-single.chip.GSEA(
									exprs = randomexprs,
									species = species[chip],
									gsdb = geneset,
									gene.set.selection = "ALL",
									sample.norm.type = "rank",
									output.score.type = "ES"
									)
    colnames(SCG)<-1:permutations
    
    save(SCG, file = paste(chip, "backgroundPerm.RData", sep = "."))
    #SCG.list<-list(SCG = SCG)
    #names(SCG.list)<-chip
    #background.SCG <- append(background.SCG, SCG.list)
  }

# this takes approx 7 days to run on the server!

# now need to produce the correlation matrix for each of these arrays

# Problem - the correlation matrix will vary between the different chip types.
# This means having to create correlation matrices for every possible chip combination!

# how much do they vary - perhaps could use an average for now.

for (i in 1:length(chipframe))
  {
   # load SCG background
   chip<-names(chipframe[i])
   print(paste("Processing", chip))
   load(paste(chip, "backgroundPerm.RData", sep = "."))
   SCG.cor<-cor(t(SCG))
   assign(paste(chip, "SCG.cor", sep = "."), SCG.cor)
  }

# plot out - N.B. pdf file too large

for (i in 1:length(chipframe))
{
  chip<-names(chipframe[i])
  print(paste("Processing", chip))
  png(paste(chip, "correlation.heatmap.png", sep = ".")) 
  heatmap(get(paste(chip, "SCG.cor", sep = ".")), scale = "none", main = chip)
  dev.off()
  }

# Assess the variation in the correlation matrix

correlation.array<-array(0, dim = c(491, 491, length(chipframe)))
for (i in 1:length(chipframe))
{
  chip<-names(chipframe[i])
  correlation.array[,,i]<-get(paste(chip, "SCG.cor", sep = "."))
  }

correlation.array.mean<-apply(correlation.array, c(1,2), mean)
correlation.array.sd<-apply(correlation.array, c(1,2), sd)
correlation.array.cv<-correlation.array.sd/abs(correlation.array.mean)

# plot cv histogram
pdf("correlation.cv.histogram.pdf")
hist(correlation.array.cv, nclass = 100, main = "Pathway Correlation Matrix", xlab = "Coefficeient of Variation")
dev.off()
# plot sd histogram
png("correlation.sd.histogram.png")
hist(correlation.array.sd, nclass = 100, main = "Pathway Correlation Matrix", xlab = "Standard deviation")
dev.off()
# plot mean against sd
png("correlation.sd.plot.png")
plot(unlist(correlation.array.mean), unlist(correlation.array.sd), main = "Pathway Correlation Matrix",
      xlab = "mean", ylab = "SD")
dev.off()

# Array that displays the greatest deviation are the GPL97 and GPL340
# this corresponds to the B chips, HGU133B and MOE430B
# could try excluding these from the analysis
Bchips<-match(c("GPL97", "GPL340"), names(chipframe))

correlation.array.Bsub.mean<-apply(correlation.array[,,-Bchips], c(1,2), mean)
correlation.array.Bsub.sd<-apply(correlation.array[,,-Bchips], c(1,2), sd)
correlation.array.Bsub.cv<-correlation.array.Bsub.sd/abs(correlation.array.Bsub.mean)
# plot mean against sd
png("correlation.Bsub.sd.plot.png")
plot(unlist(correlation.array.Bsub.mean), unlist(correlation.array.Bsub.sd), main = "Pathway Correlation Matrix",
      xlab = "mean", ylab = "SD")
dev.off()
# plot cv histogram
pdf("correlation.Bsub.cv.histogram.pdf")
hist(correlation.array.Bsub.cv, nclass = 100, main = "Pathway Correlation Matrix", xlab = "Coefficeient of Variation")
dev.off()

# output the averaged matrix
# N.B. this didn't save first time AND WAS INCORRECT AS CALCULATED BASED ON CV NOT MEAN
#pathwayCorrelationMatrix<-correlation.array.Bsub.cv
colnames(pathwayCorrelationMatrix)<-names(kegg_wiki_TR_static.Hs.gs)
rownames(pathwayCorrelationMatrix)<-names(kegg_wiki_TR_static.Hs.gs)

write(pathwayCorrelationMatrix, file = "pathwayCorrelationMatrix.RData")

save.image("covarianceBackground.RData")

#############
# remake correlation matrix

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data")
# load genesets
load("/home/galtschu2/Documents/Databases/Gene sets/kegg_wiki_TR_static.Hs.gs")
load("/home/galtschu2/Documents/Databases/Gene sets/kegg_wiki_TR_static.Mm.gs")



# load chip data
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/chipframe_update.RData")

for (i in 1:length(chipframe))
{
 # load SCG background
 chip<-names(chipframe[i])
 print(paste("Processing", chip))
 load(paste(chip, "backgroundPerm.RData", sep = "."))
 SCG.cor<-cor(t(SCG))
 assign(paste(chip, "SCG.cor", sep = "."), SCG.cor)
}


correlation.array<-array(0, dim = c(491, 491, length(chipframe)))
for (i in 1:length(chipframe))
{
  chip<-names(chipframe[i])
  correlation.array[,,i]<-get(paste(chip, "SCG.cor", sep = "."))
  }

correlation.array.mean<-apply(correlation.array, c(1,2), mean)
correlation.array.sd<-apply(correlation.array, c(1,2), sd)
correlation.array.cv<-correlation.array.sd/abs(correlation.array.mean)

Bchips<-match(c("GPL97", "GPL340"), names(chipframe))

correlation.array.Bsub.mean<-apply(correlation.array[,,-Bchips], c(1,2), mean)
correlation.array.Bsub.sd<-apply(correlation.array[,,-Bchips], c(1,2), sd)
correlation.array.Bsub.cv<-correlation.array.Bsub.sd/abs(correlation.array.Bsub.mean)


# output the averaged matrix

pathwayCorrelationMatrix<-correlation.array.Bsub.mean
colnames(pathwayCorrelationMatrix)<-names(kegg_wiki_TR_static.Hs.gs)
rownames(pathwayCorrelationMatrix)<-names(kegg_wiki_TR_static.Hs.gs)

save(pathwayCorrelationMatrix, file = "pathwayCorrelationMatrix_2.RData")
