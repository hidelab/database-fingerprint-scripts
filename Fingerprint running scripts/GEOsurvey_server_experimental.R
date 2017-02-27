# Updating fingerprint database with exprimental datasets
# Author: Gabriel Altschuler
# Status: Refactoring
# Timestamp: 20110413
# Run on hpc111 by source("/home/galtschu2/database-fingerprint-scripts/Fingerprint running scripts/GEOsurvey_server_experimental.R")
# This script processes and saves fingeprints for the experimental genelists

# define current directory
running.dir <- dirname(parent.frame(2)$ofile)
base.scripts<- gsub("database-fingerprint-scripts/Fingerprint running scripts", "fingerprint/scripts", running.dir)
dataPath<- gsub("database-fingerprint-scripts/Fingerprint running scripts", "fingerprint/data", running.dir)
if(running.dir == "."){
	base.scripts <- "../../fingerprint/scripts"
	dataPath <- "../../fingerprint/data"
	}
	
source(paste(base.scripts, "geo2fingerprint.R", sep="/"))

if (file.exists("Reading in samples_server.R")) {
  source('Reading in samples_server.R')  
} else {
  source(paste(running.dir, "Reading in samples_server.R", sep="/"))
  }


if (!(exists("chipframe"))){
	load(paste(dataPath, "chipframe.RData", sep="/"), .GlobalEnv)
	}
if (!(exists("ExperimentalGeneSigs.Hs.gs"))){
	load(paste(dataPath, "ExperimentalGeneSigs.Hs.gs", sep="/"), .GlobalEnv)
	}
if (!(exists("ExperimentalGeneSigs.Mm.gs"))){
	load(paste(dataPath, "ExperimentalGeneSigs.Mm.gs", sep="/"), .GlobalEnv)
	}
setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
library(foreach)
library(doMC)
registerDoMC(cores = 14)



# run scripts to process using the experimental options
a=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_a_experimental.R", sep="/")))
b=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_b_experimental.R", sep="/")))
c=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_c_experimental.R", sep="/")))
d=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_d_experimental.R", sep="/")))
e=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_e_experimental.R", sep="/")))
f=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_f_experimental.R", sep="/")))
g=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_g_experimental.R", sep="/")))
h=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_h_experimental.R", sep="/")))
i=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_i_experimental.R", sep="/")))
j=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_j_experimental.R", sep="/")))
k=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_k_experimental.R", sep="/")))
l=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_l_experimental.R", sep="/")))
m=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_m_experimental.R", sep="/")))
n=parallel(source(paste(running.dir, "Experimental_parallel", "parallel_n_experimental.R", sep="/")))

		
# Collect processes at the end to terminate jobs
Sys.sleep(10)
x=collect(list(a,b,c,d,e,f,g,h,i,j,k,l,m,n))


