# produce random genesets to test fingerprint
# Author: Gabriel Altschuler
# Timestamp: 20110611
# Status: Completed
# Run location: hpc111
# Addressing Lincoln Stein's question
##
# ""Do fingerprints based on the FI network modules work
# better than arbitrarily dividing the genome
# into a like number and distribution of random modules.
# That is, is partitioning by pathway information actually
# contributing to the accuracy of the fingerprint matching
# or would other partitionings work as well?""
##

# Need to create an arbitrary set of pathways
# These should be
# 1) based on the same genes
# 2) have the same distribution of pathway sizes
# N.B. the correlation structure will be lost
# This would be too complicated to account for

# Easiest way to achieve this is to 
# a) turn all pathways into one long vector
# b) shuffle the order
# c) Re-partition the list into the same pathways


source("/home/galtschu2/fingerprint/scripts/gabriel functions.R")

# load geneset
load("/home/galtschu2/fingerprint/data/pathprint.v0.3.Hs.gs")

gene.vector<-unlist(pathprint.v0.3.Hs.gs)
gene.vector<-sample(gene.vector)
random.v0.3.Hs.gs<-relist(flesh = gene.vector, skeleton = pathprint.v0.3.Hs.gs)
names(random.v0.3.Hs.gs)<-paste("PERMUTED", names(random.v0.3.Hs.gs), sep = "_")
 
# now convert for each of the additional species
homologene.dir<-"/home/galtschu2/Documents/Databases/Homologene"
setwd(homologene.dir)
taxID<-read.delim("taxid_taxname", header = FALSE, stringsAsFactors = FALSE)
homologeneFile<-"/home/galtschu2/Documents/Databases/Homologene/homologene.data"
# list of additional species (N.B. no pig)
additional.species<-c("Mus musculus", "Rattus norvegicus", "Drosophila melanogaster", "Danio rerio", "Caenorhabditis elegans")

initialize<-function(x){
     part.2<-unlist(strsplit(x, " "))[2]
     initials<-paste(substring(x, 1, 1),
   	substring(part.2, 1, 1),
 		sep = "")
     return(initials)
     }

additional.species.initials<-sapply(additional.species, initialize)

# create genesets for other species
for (i in 1:length(additional.species)){
  geneset<-genesetConvert(random.v0.3.Hs.gs, convertTo = taxID[match(additional.species[i], taxID[,2]), 1], homologeneFile = homologeneFile)
  assign(paste("random.v0.3", additional.species.initials[i], "gs", sep = "."), geneset)
  }

# save
save(random.v0.3.Hs.gs, file = "/home/galtschu2/fingerprint/data/random.v0.3.Hs.gs")
for (i in 1:length(additional.species)){
  geneset<-paste("random.v0.3", additional.species.initials[i], "gs", sep = ".")
  save(list = geneset, file = paste("/home/galtschu2/fingerprint/data/", geneset, sep = ""))
  }

# Additional part, run locally to create look-up for which geneset to use

species<-c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Drosophila melanogaster", "Danio rerio", "Caenorhabditis elegans")
initialize<-function(x){
     part.2<-unlist(strsplit(x, " "))[2]
     initials<-paste(substring(x, 1, 1),
     substring(part.2, 1, 1),
 		sep = "")
     return(initials)
     }
species.initials<-sapply(species, initialize)

genesets<-sapply(species.initials, function(x){paste("random.v0.3", x, "gs", sep = ".")})

# quick fix to add common names
genesets<-c(genesets, genesets)
names(genesets)<-c(names(genesets)[1:6], c("human", "mouse", "rat", "drosophila", "zebrafish", "C.elegans"))

save(genesets, file = "/Users/GabrielAltschuler/Dropbox/fingerprint/data/random.genesets.RData")



load("/Users/GabrielAltschuler/Dropbox/fingerprint/data/random.genesets.RData")


# end of script




# end of script

