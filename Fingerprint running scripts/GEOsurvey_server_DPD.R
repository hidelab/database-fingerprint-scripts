# Running fingerprint based on squared rank enrichment scores using Drug, pathway, disease genesets compiled by Gabriel
# Author: Gabriel Altschuler
# Timestamp: 20121117
# Status: Construction
# Based on pathprint package 1.2.2
# Loads data from the pathprint dataframe, not GEO
setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
# load packages
library(GMAfunctions)
library(pathprint)

# pathprint lazyloads the chipframe (chip annotation details) and genesets.

# need to overwrite pathprint genesets with the genesets that correspond to the modules
setwd("/home/galtschu2/database-fingerprint-scripts/DPDgenesets")
load("DPD.genesets.RData")
lapply(genesets, load, envir = .GlobalEnv)
setwd("/home/galtschu2/Documents/Projects/Fingerprinting")

# source a script to read the list of GEO files to be fingerprinted
# This creates an object, sampleset, which is a list containing the full list of GSMs
# The list is randomly split into 14 parts, to be distributed across the processors
# only select human HGU133PLUS2 data
source("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/Reading in samples_NO_UPDATE_GPL570.R")

library(foreach)
library(doMC)
registerDoMC(cores = 15)

# define path directory, and create if necessary
path<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/DPD/"
try(system(paste("mkdir", path, sep = " ")))

GEOpath<-"/home/galtschu2/Documents/Databases/GEOfiles/"

setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts")

# these scripts run for single chip enrichment - mean rank not mean squared rank
a=parallel(source("DPD_parallel/DPD_parallel_a.R"))
b=parallel(source("DPD_parallel/DPD_parallel_b.R"))
c=parallel(source("DPD_parallel/DPD_parallel_c.R"))
d=parallel(source("DPD_parallel/DPD_parallel_d.R"))
e=parallel(source("DPD_parallel/DPD_parallel_e.R"))
f=parallel(source("DPD_parallel/DPD_parallel_f.R"))
g=parallel(source("DPD_parallel/DPD_parallel_g.R"))
h=parallel(source("DPD_parallel/DPD_parallel_h.R"))
i=parallel(source("DPD_parallel/DPD_parallel_i.R"))
j=parallel(source("DPD_parallel/DPD_parallel_j.R"))
k=parallel(source("DPD_parallel/DPD_parallel_k.R"))
l=parallel(source("DPD_parallel/DPD_parallel_l.R"))
m=parallel(source("DPD_parallel/DPD_parallel_m.R"))
n=parallel(source("DPD_parallel/DPD_parallel_n.R"))

# collect processes
Sys.sleep(10)
x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))



