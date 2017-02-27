# Conversion of gene sets to multiple species
# Author: Gabriel Altschuler
# Timestamp: 06082011
# Status: Completed

# download homologene file to server
# ftp ftp.ncbi.nih.gov
# cd /pub/HomoloGene/current
# get homologene.data
# cd build_inputs
# get taxid_taxname
# bye

source("/home/galtschu2/fingerprint/scripts/gabriel functions.R")

# load current datasource

library(pathprint)
data(kegg_wiki_TR_static.Hs.gs)
pathprint.v0.2.Hs.gs<-kegg_wiki_TR_static.Hs.gs

# also load next geneset list
load("/home/galtschu2/fingerprint/data/reactome.kegg.wiki.netpath.static.Hs.gs")
pathprint.v1.Hs.gs<-reactome.kegg.wiki.netpath.static.Hs.gs

# load combined geneset list
load("/home/galtschu2/fingerprint/data/pathprint.v0.3.Hs.gs")

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

for (i in 1:length(additional.species)){
  geneset<-genesetConvert(pathprint.v0.2.Hs.gs, convertTo = taxID[match(additional.species[i], taxID[,2]), 1], homologeneFile = homologeneFile)
  assign(paste("pathprint.v0.2", additional.species.initials[i], "gs", sep = "."), geneset)
  }

# save
save(pathprint.v0.2.Hs.gs, file = "/home/galtschu2/fingerprint/data/pathprint.v0.2.Hs.gs")
for (i in 1:length(additional.species)){
  geneset<-paste("pathprint.v0.2", additional.species.initials[i], "gs", sep = ".")
  save(list = geneset, file = paste("/home/galtschu2/fingerprint/data/", geneset, sep = ""))
  }


# pathprint v.1 list
for (i in 1:length(additional.species)){
  geneset<-genesetConvert(pathprint.v1.Hs.gs, convertTo = taxID[match(additional.species[i], taxID[,2]), 1], homologeneFile = homologeneFile)
  assign(paste("pathprint.v1", additional.species.initials[i], "gs", sep = "."), geneset)
  }

# save
save(pathprint.v1.Hs.gs, file = "/home/galtschu2/fingerprint/data/pathprint.v1.Hs.gs")
for (i in 1:length(additional.species)){
  geneset<-paste("pathprint.v1", additional.species.initials[i], "gs", sep = ".")
  save(list = geneset, file = paste("/home/galtschu2/fingerprint/data/", geneset, sep = ""))
  }


# pathprint v.03 list
for (i in 1:length(additional.species)){
  geneset<-genesetConvert(pathprint.v0.3.Hs.gs, convertTo = taxID[match(additional.species[i], taxID[,2]), 1], homologeneFile = homologeneFile)
  assign(paste("pathprint.v0.3", additional.species.initials[i], "gs", sep = "."), geneset)
  }

# save
for (i in 1:length(additional.species)){
  geneset<-paste("pathprint.v0.3", additional.species.initials[i], "gs", sep = ".")
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

genesets<-sapply(species.initials, function(x){paste("pathprint.v0.3", x, "gs", sep = ".")})

# quick fix to add common names
genesets<-c(genesets, genesets)
names(genesets)<-c(names(genesets)[1:6], c("human", "mouse", "rat", "drosophila", "zebrafish", "C.elegans"))

save(genesets, file = "/Users/GabrielAltschuler/Dropbox/fingerprint/data/genesets.RData")



load("/Users/GabrielAltschuler/Dropbox/fingerprint/data/genesets.RData")


# end of script
