# Produce fingerprint scores based on random permutations_using v0.3 pathways
# Author: Gabriel Altschuler
# Status: Updated for 52057 permutations
# Timestamp: 20110726

# Attempt produce the SCE score based on 52057 randomly permuted genesets
# Just perform for GPL570

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data")

#######
# Edit for increased number of permutations
# load("GPL570.GEO.v0.3.RData")
# ncol(GPL570.GEO.v0.3)
# 52057

# already have 10000 permutations in the previous dataframe

load("/home/galtschu2/fingerprint/data/pathprint.v0.3.Hs.gs")
load("/home/galtschu2/fingerprint/data/chipframe.RData")

# need to perform 42057 permutations. This takes up too much memory.
# This is in-part due to inefficiencies in the single.chip.GSEA routine
# Instead run in 9 chucks of 4673
# Should really parallelize this but if not a priority then not needed

permutations <- 4673
chipdata<-chipframe[["GPL570"]]$ann
chipgenes<-unique(sub("_at", "", chipdata$EntrezID))

# produce random array, dim (Ngenes x Npermutations)
for (i in 1:9){
  print(paste("part", i, "of 9", sep = " "))

  randomexprs<-apply(matrix(rep(1:length(chipgenes),permutations),
                            ncol = permutations
                             ),
                     2, sample
                     )
  rownames(randomexprs)<-chipgenes
  library(pathprint)
  GPL570.random.v0.3<-single.chip.enrichment(
    				exprs = randomexprs,
  					geneset = pathprint.v0.3.Hs.gs,
  					transformation = "rank",
  					statistic = "mean",
  					normalizedScore = FALSE,
  					progressBar = TRUE
  					)
  #save(GPL570.random.v0.3, file = "GPL570.random.v0.3.RData")
  save(GPL570.random.v0.3, file = paste("GPL570.random.", i, ".v0.3.RData"))
}

# now re-assemble into a large dataframe
GPL570.random.v0.3.52057<-matrix(nrow = nrow(GPL570.random.v0.3), ncol = 52057)
rownames(GPL570.random.v0.3.52057)<-rownames(GPL570.random.v0.3)

for (i in 1:9){
  load(paste("GPL570.random.", i, ".v0.3.RData"))
  start<-1+((i-1)*4673)
  end <- i*4673
  GPL570.random.v0.3.52057[,start:end]<-GPL570.random.v0.3
}
load("GPL570.random.v0.3.RData")
GPL570.random.v0.3.52057[,42058:52057]<-GPL570.random.v0.3
save(GPL570.random.v0.3.52057, file = "GPL570.random.v0.3.52057.RData")
# Comparison is with the SCE matrix for GPL570 
# load files

load("/home/galtschu2/Documents/Projects/Fingerprinting/data/platform.frame.2011-06-12.RData")

system("ls /home/galtschu2/Documents/Projects/Fingerprinting/data/v0.3_.frame*")
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/v0.3_.frame.2011-06-12.RData")

platform.frame.GPL570<-platform.frame[platform.frame$Platform == "GPL570",]

GPL570.GEO.v0.3<-SCE.frame[match(platform.frame.GPL570$GEO, rownames(SCE.frame)),]
GPL570.GEO.v0.3<-t(GPL570.GEO.v0.3)

save(GPL570.GEO.v0.3, file = "GPL570.GEO.v0.3.RData")

# can't threshold as don't have the threshold values yet
#GPL570.GEO.fingerprint<-thresholdFingerprint(GPL570.GEO.v0.3, "GPL570")
#save(GPL570.GEO.fingerprint, file = "GPL570.GEO.fingerprint.RData")