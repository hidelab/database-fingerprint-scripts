# Running GEO fingerprinting using random genesets

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

# LOADING RANDOM GENESETS
load("/home/galtschu2/fingerprint/data/random.genesets.RData")
sapply(genesets, function(x){load(paste("/home/galtschu2/fingerprint/data/", x, sep = ""), .GlobalEnv)})
print("Running random genesets")
print(genesets)

	
library(foreach)
library(doMC)
registerDoMC(cores = 15)


# these scripts run for single chip enrichment using random genesets
a=parallel(source("random_SCE_parallel/random_SCE_parallel_a.R"))
b=parallel(source("random_SCE_parallel/random_SCE_parallel_b.R"))
c=parallel(source("random_SCE_parallel/random_SCE_parallel_c.R"))
d=parallel(source("random_SCE_parallel/random_SCE_parallel_d.R"))
e=parallel(source("random_SCE_parallel/random_SCE_parallel_e.R"))
f=parallel(source("random_SCE_parallel/random_SCE_parallel_f.R"))
g=parallel(source("random_SCE_parallel/random_SCE_parallel_g.R"))
h=parallel(source("random_SCE_parallel/random_SCE_parallel_h.R"))
i=parallel(source("random_SCE_parallel/random_SCE_parallel_i.R"))
j=parallel(source("random_SCE_parallel/random_SCE_parallel_j.R"))
k=parallel(source("random_SCE_parallel/random_SCE_parallel_k.R"))
l=parallel(source("random_SCE_parallel/random_SCE_parallel_l.R"))
m=parallel(source("random_SCE_parallel/random_SCE_parallel_m.R"))
n=parallel(source("random_SCE_parallel/random_SCE_parallel_n.R"))

# collect processes
Sys.sleep(10)
x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))


