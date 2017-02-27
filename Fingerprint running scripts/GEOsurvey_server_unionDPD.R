# Running fingerprint based on median rank enrichment scores using DPD genesets compiled by Gabriel
# Uses genesets with consistent representation between mouse GPL1261 and human GPL570
# Author: Gabriel Altschuler
# Timestamp: 20121117
# Status: Construction
# Loads data from the pathprint dataframe, not GEO
setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
# load packages
library(GMAfunctions)
library(pathprint)

# pathprint lazyloads the chipframe (chip annotation details) and genesets.

# need to overwrite pathprint genesets with the genesets that correspond to the modules
setwd("/home/galtschu2/database-fingerprint-scripts/unionGenesets")
load("DPDUnion.genesets.RData")
lapply(genesets, load, envir = .GlobalEnv)
# load chipframe with alternative GPL570 and GPL1261 annotations using only the commom genes
load("chipframeCommon_GPL570_GPL1261.RData")
setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
header = "UnionDPD"
testStat = "median"
# source a script to read the list of GEO files to be fingerprinted
# This creates an object, sampleset, which is a list containing the full list of GSMs
# The list is randomly split into 14 parts, to be distributed across the processors
# only select human HGU133PLUS2 data
source("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/Reading in samples_NO_UPDATE_GPL570_GPL1261.R")

library(foreach)
library(doMC)
registerDoMC(cores = 15)

# define path directory, and create if necessary
path<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/UnionDPD/"
try(system(paste("mkdir", path, sep = " ")))

GEOpath<-"/home/galtschu2/Documents/Databases/GEOfiles/"

setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts")

# these scripts run for single chip enrichment - median rank not mean squared rank
a=parallel(source("parallel/parallel_a.R"))
b=parallel(source("parallel/parallel_b.R"))
c=parallel(source("parallel/parallel_c.R"))
d=parallel(source("parallel/parallel_d.R"))
e=parallel(source("parallel/parallel_e.R"))
f=parallel(source("parallel/parallel_f.R"))
g=parallel(source("parallel/parallel_g.R"))
h=parallel(source("parallel/parallel_h.R"))
i=parallel(source("parallel/parallel_i.R"))
j=parallel(source("parallel/parallel_j.R"))
k=parallel(source("parallel/parallel_k.R"))
l=parallel(source("parallel/parallel_l.R"))
m=parallel(source("parallel/parallel_m.R"))
n=parallel(source("parallel/parallel_n.R"))

# collect processes
Sys.sleep(10)
x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))



