setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
source("/home/galtschu2/fingerprint/scripts/gabriel functions.R")

# need to load pathprint for annotation script but load updated chipframe on top
library(pathprint)
load("/home/galtschu2/fingerprint/data/chipframe.RData", .GlobalEnv)

# source script that reads the list of GEO files to be fingerprinted
# This creates an object, sampleset, which is a list containing the full list of GSMs
# The list is randomly split into 14 parts, to be distributed across the processors
source("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/Reading in samples_server.R")

# Now load all the dataframes required for the processing of the data from each chip
# Load into the global environment so available for all processes
setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts")
#if (!(exists("chipframe"))){
#	load("/home/galtschu2/fingerprint/data/chipframe.RData", .GlobalEnv)
#	}
load("/home/galtschu2/fingerprint/data/genesets.RData")

sapply(genesets, function(x){load(paste("/home/galtschu2/fingerprint/data/", x, sep = ""), .GlobalEnv)})

#if (!(exists("kegg_wiki_TR_static.Hs.gs"))){
#	load("/home/galtschu2/Fingerprint/data/kegg_wiki_TR_static.Hs.gs", .GlobalEnv)
#	}
#if (!(exists("kegg_wiki_TR_static.Mm.gs"))){
#	load("/home/galtschu2/Fingerprint/data/kegg_wiki_TR_static.Mm.gs", .GlobalEnv)
#	}
	
	
library(foreach)
library(doMC)
registerDoMC(cores = 15)

# these ran for the single chip GSEA
# a=parallel(source("scripts/parallel_a.R"))
# b=parallel(source("scripts/parallel_b.R"))
# c=parallel(source("scripts/parallel_c.R"))
# d=parallel(source("scripts/parallel_d.R"))
# e=parallel(source("scripts/parallel_e.R"))
# f=parallel(source("scripts/parallel_f.R"))
# g=parallel(source("scripts/parallel_g.R"))
# h=parallel(source("scripts/parallel_h.R"))
# i=parallel(source("scripts/parallel_i.R"))
# j=parallel(source("scripts/parallel_j.R"))
# k=parallel(source("scripts/parallel_k.R"))
# l=parallel(source("scripts/parallel_l.R"))
# m=parallel(source("scripts/parallel_m.R"))
# n=parallel(source("scripts/parallel_n.R"))

# these scripts run for single chip enrichment
a=parallel(source("SCE_parallel/SCE_parallel_a.R"))
b=parallel(source("SCE_parallel/SCE_parallel_b.R"))
c=parallel(source("SCE_parallel/SCE_parallel_c.R"))
d=parallel(source("SCE_parallel/SCE_parallel_d.R"))
e=parallel(source("SCE_parallel/SCE_parallel_e.R"))
f=parallel(source("SCE_parallel/SCE_parallel_f.R"))
g=parallel(source("SCE_parallel/SCE_parallel_g.R"))
h=parallel(source("SCE_parallel/SCE_parallel_h.R"))
i=parallel(source("SCE_parallel/SCE_parallel_i.R"))
j=parallel(source("SCE_parallel/SCE_parallel_j.R"))
k=parallel(source("SCE_parallel/SCE_parallel_k.R"))
l=parallel(source("SCE_parallel/SCE_parallel_l.R"))
m=parallel(source("SCE_parallel/SCE_parallel_m.R"))
n=parallel(source("SCE_parallel/SCE_parallel_n.R"))

# collect processes
Sys.sleep(10)
#x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l))
x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))


