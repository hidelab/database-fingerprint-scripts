# Fixing bugs in running script

# halt scripts using ^C

# see which have crahsed early
processes<-letters[1:14]
for(i in 1:length(processes)){
  collect(get(processes[i]))
}

collect(list(a,b,c,d,e,f,h,i,j,k,l,m,n))

setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/v0.3/")

fingerprinted<-gsub("v0.3_", "", dir())
# find first stop
GEOsamples.processed<-GEOsamples.rand %in% fingerprinted
names(GEOsamples.processed)<-GEOsamples.rand
GEOsamples.processed[match(FALSE, GEOsamples.processed)]


head(GEOsamples.processed[GEOsamples.processed == FALSE])


# test arrays

temp<-geo2fingerprint(
              GSM = "GSM79476", 
              GEOthreshold = FALSE,
              GEOpath = "/home/galtschu2/Documents/Databases/GEOfiles/",
              geneset = "KEGG and Wikipathways and static",
              enrichmentMethod = "SCE",
              transformation = "rank",
              statistic = "mean",
              normalizedScore = FALSE,
              progressBar = FALSE
              )
              
# fix bugs

# now clean up
#kill processes if not done by collect already

processes<-letters[1:14]
process.id<-vector("numeric", 14)
for(i in 1:length(processes)){
  process.id<-names(collect(get(processes[i])))
  system(paste("kill -9", process.id, sep = " "))
  }

# now quit and restart fingerprinting 

# post-processing bug - found that there are some chips that give NAs due to incomplete download of GEO files
# need to find these and re-fingerprint
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")
system("ls v0.3_.frame*")
load("v0.3_.frame.2011-06-12.RData")
na<-which(is.na(SCE.frame), arr.ind = TRUE)
na.arrays<-unique(rownames(na))
# there are two possible sources for these NAs
# either there is no representation for a pathway on that particular chip
# or the GEO record was not properly downloaded for that chip
# check to see which platforms they come from
system("ls platform.frame*")
load("platform.frame.2011-06-12.RData")

na.platforms<-platform.frame[na.arrays,]
table(na.platforms$Platform)
table(platform.frame$Platform)

# compare the size of a typical vs an atypical GEO record
na.platforms$GEO[match("GPL1261", na.platforms$Platform)]
platform.frame$GEO[match("GPL1261", platform.frame$Platform)]

na.platforms$GEO[match("GPL570", na.platforms$Platform)]
platform.frame$GEO[match("GPL570", platform.frame$Platform)]
head(na.platforms$GEO[na.platforms$Platform = "GPL570"])

system("ls /home/galtschu2/Documents/Databases/GEOfiles/GSM303303* -l")
system("ls /home/galtschu2/Documents/Databases/GEOfiles/GSM100169* -l")
system("ls /home/galtschu2/Documents/Databases/GEOfiles/GSM132948* -l")
system("ls /home/galtschu2/Documents/Databases/GEOfiles/GSM100200* -l")

# compare local version with newly downloaded
GSM132948.local<-getGEO("GSM132948", destdir = "/home/galtschu2/Documents/Databases/GEOfiles/")
GSM132948<-getGEO("GSM132948")

GSM303303.local<-getGEO("GSM303303", destdir = "/home/galtschu2/Documents/Databases/GEOfiles/")
GSM303303 <-getGEO("GSM303303")
GSM100200.local<-getGEO("GSM100200", destdir = "/home/galtschu2/Documents/Databases/GEOfiles/")
GSM100200 <-getGEO("GSM100200")
GSM400264.local<-getGEO("GSM400264", destdir = "/home/galtschu2/Documents/Databases/GEOfiles/")
GSM400264 <-getGEO("GSM400264")


# There is a disparity
# The most direct way to remedy this is to remove all of these files and re-fingerprint
# but this is not always the reason - try to track this down

library(pathprint)
data(kegg_wiki_TR_static.Hs.gs)
data(chipframe)
GSM132948.SCE<-single.chip.enrichment(customCDFAnn(Table(GSM132948), chipframe$GPL570$ann), geneset = kegg_wiki_TR_static.Hs.gs)
GSM132948.local.SCE<-single.chip.enrichment(customCDFAnn(Table(GSM132948.local), chipframe$GPL570$ann), geneset = kegg_wiki_TR_static.Hs.gs)

GSM400264.SCE<-single.chip.enrichment(customCDFAnn(Table(GSM400264), chipframe$GPL570$ann), geneset = kegg_wiki_TR_static.Hs.gs)
GSM400264.local.SCE<-single.chip.enrichment(customCDFAnn(Table(GSM400264.local), chipframe$GPL570$ann), geneset = kegg_wiki_TR_static.Hs.gs)

head(SCE.frame["GSM400264",])

load("/home/galtschu2/fingerprint/data/genesets.RData")
sapply(genesets, function(x){load(paste("/home/galtschu2/fingerprint/data/", x, sep = ""), .GlobalEnv)})

temp<-geo2fingerprint(
              GSM = "GSM400264", 
              GEOthreshold = FALSE,
              GEOpath = "/home/galtschu2/Documents/Databases/GEOfiles/",
              geneset = "KEGG and Wikipathways and static",
              enrichmentMethod = "SCE",
              transformation = "rank",
              statistic = "mean",
              normalizedScore = FALSE,
              progressBar = FALSE
              )
exprs<-GSMtable2exprs(Table(GSM400264))
geo.SCG <- exprs2fingerprint_options(											exprs = exprs, platform = "GPL570", species = "human", GEOthreshold = FALSE, geneset = pathprint.v0.3.Hs.gs,normalizedScore = FALSE
										)	

# this one is due to a human/mouse discrepancy in the metadata