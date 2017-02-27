# Produce fingerprint scores based on random permutations
# Author: Gabriel Altschuler
# Status: Completed
# Timestamp: 20110609

# Attempt produce the SCE score based on 10000 randomly permuted genesets
# Just perform for GPL570

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data")

library(pathprint)
data(kegg_wiki_TR_static.Hs.gs)
data(chipframe)

permutations <- 10000
chipdata<-chipframe[["GPL570"]]$ann
chipgenes<-unique(sub("_at", "", chipdata$EntrezID))

# produce random array, dim (Ngenes x Npermutations)
randomexprs<-apply(matrix(rep(1:length(chipgenes),permutations),
                          ncol = permutations
                           ),
                   2, sample
                   )
rownames(randomexprs)<-chipgenes

GPL570.random.SCE<-single.chip.enrichment(
					exprs = randomexprs,
					geneset = kegg_wiki_TR_static.Hs.gs,
					transformation = "rank",
					statistic = "mean",
					normalizedScore = FALSE,
					progressBar = TRUE
					)
save(GPL570.random.SCE, file = "GPL570.random.SCE.RData")

# produce a ternary fingerprint based on these values

GPL570.random.fingerprint<-thresholdFingerprint(GPL570.random.SCE, "GPL570")
save(GPL570.random.fingerprint, file = "GPL570.random.fingerprint.RData")

# don't think you can really compare the finerprint values as this is built on the biological null distributions rather than the random distributions
# better comparison is with the SCE matrix for GPL570 
# load files

load("/home/galtschu2/Documents/Projects/Fingerprinting/data/pathwayFrames/platform.frame.2011-04-05.RData")

system("ls /home/galtschu2/Documents/Projects/Fingerprinting/data/SCG.frame*")
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/SCG.frame.2011-04-12.RData")

platform.frame.GPL570<-platform.frame[platform.frame$Platform == "GPL570",]

GPL570.GEO.SCE<-SCG.frame[match(platform.frame.GPL570$GEO, rownames(SCG.frame)),]
GPL570.GEO.SCE<-t(GPL570.GEO.SCE)
colnames(GPL570.GEO.SCE)<-gsub("SCE_", "", colnames(GPL570.GEO.SCE))
colnames(GPL570.GEO.SCE)<-gsub("SCE_", "", colnames(GPL570.GEO.SCE))

save(GPL570.GEO.SCE, file = "GPL570.GEO.SCE.RData")
GPL570.GEO.fingerprint<-thresholdFingerprint(GPL570.GEO.SCE, "GPL570")
save(GPL570.GEO.fingerprint, file = "GPL570.GEO.fingerprint.RData")