# Running fingerprint based on squared-rank enrichment scores
# Author: Gabriel Altschuler
# Timestamp: 11282011
# Status: Updated
# Based on pathprint package 1.2

setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
# load packages
library(GMAfunctions)
library(pathprint)

# pathprint lazyloads the chipframe (chip annotation details) and genesets.

# source a script to read the list of GEO files to be fingerprinted
# This creates an object, sampleset, which is a list containing the full list of GSMs
# The list is randomly split into 14 parts, to be distributed across the processors
source("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/Reading in samples_server.R")

library(foreach)
library(doMC)
registerDoMC(cores = 15)

# define path directory
path<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints/sq/"
GEOpath<-"/home/galtschu2/Documents/Databases/GEOfiles/"

setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts")

# these scripts run for single chip enrichment
a=parallel(source("SCE_sq_parallel/SCE_parallel_a.R"))
b=parallel(source("SCE_sq_parallel/SCE_parallel_b.R"))
c=parallel(source("SCE_sq_parallel/SCE_parallel_c.R"))
d=parallel(source("SCE_sq_parallel/SCE_parallel_d.R"))
e=parallel(source("SCE_sq_parallel/SCE_parallel_e.R"))
f=parallel(source("SCE_sq_parallel/SCE_parallel_f.R"))
g=parallel(source("SCE_sq_parallel/SCE_parallel_g.R"))
h=parallel(source("SCE_sq_parallel/SCE_parallel_h.R"))
i=parallel(source("SCE_sq_parallel/SCE_parallel_i.R"))
j=parallel(source("SCE_sq_parallel/SCE_parallel_j.R"))
k=parallel(source("SCE_sq_parallel/SCE_parallel_k.R"))
l=parallel(source("SCE_sq_parallel/SCE_parallel_l.R"))
m=parallel(source("SCE_sq_parallel/SCE_parallel_m.R"))
n=parallel(source("SCE_sq_parallel/SCE_parallel_n.R"))

# collect processes
Sys.sleep(10)
x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))




#################################################
## Below is the previous version of the script ##
#################################################

#setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
##source("/home/galtschu2/fingerprint/scripts/gabriel functions.R")
#
#
## need to load pathprint for annotation script but load updated chipframe on top
#library(pathprint.v0.3.beta3)
#data(chipframe)
#
## source script that reads the list of GEO files to be fingerprinted
## This creates an object, sampleset, which is a list containing the full list of GSMs
## The list is randomly split into 14 parts, to be distributed across the processors
#source("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/Reading in samples_server.R")
#
## Now load all the dataframes required for the processing of the data from each chip
## Load into the global environment so available for all processes
#setwd("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts")
##if (!(exists("chipframe"))){
##	load("/home/galtschu2/fingerprint/data/chipframe.RData", .GlobalEnv)
##	}
#load("/home/galtschu2/fingerprint/data/genesets.RData")
#
#sapply(genesets, function(x){load(paste("/home/galtschu2/fingerprint/data/", x, sep = ""), .GlobalEnv)})
#
#
#	
#library(foreach)
#library(doMC)
#registerDoMC(cores = 15)
#
#
## these scripts run for single chip enrichment
#a=parallel(source("SCE_sq_parallel/SCE_parallel_a.R"))
#b=parallel(source("SCE_sq_parallel/SCE_parallel_b.R"))
#c=parallel(source("SCE_sq_parallel/SCE_parallel_c.R"))
#d=parallel(source("SCE_sq_parallel/SCE_parallel_d.R"))
#e=parallel(source("SCE_sq_parallel/SCE_parallel_e.R"))
#f=parallel(source("SCE_sq_parallel/SCE_parallel_f.R"))
#g=parallel(source("SCE_sq_parallel/SCE_parallel_g.R"))
#h=parallel(source("SCE_sq_parallel/SCE_parallel_h.R"))
#i=parallel(source("SCE_sq_parallel/SCE_parallel_i.R"))
#j=parallel(source("SCE_sq_parallel/SCE_parallel_j.R"))
#k=parallel(source("SCE_sq_parallel/SCE_parallel_k.R"))
#l=parallel(source("SCE_sq_parallel/SCE_parallel_l.R"))
#m=parallel(source("SCE_sq_parallel/SCE_parallel_m.R"))
#n=parallel(source("SCE_sq_parallel/SCE_parallel_n.R"))
#
## collect processes
#Sys.sleep(10)
##x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l))
#x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))
#

