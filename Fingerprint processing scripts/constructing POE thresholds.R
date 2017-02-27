# Constructs the POE thresholds for the fingerprint
#
# Author: Gabriel Altschuler
# Status: Updated to allow source
# Timestamp: 20111201
# Location: hpc111
# Script to obtain threshold values for the fingerprint
# These values will be passed to a dataframe, thresholds.RData
# This will be included into the fingerprinting R package
# The background parameters are determined by the PEO method
# The threshold was optimized against a panel of tissue samples

# define path of repository
# pathprintRepository<-"/home/galtschu2/fingerprint/data/"
definePath<-function(){
readline("define pathprint repository, or blank for default (/home/galtschu2/fingerprint/data/) : ")
}
pathprintRepository<-definePath()
if (pathprintRepository == ""){
	pathprintRepository<-"/home/galtschu2/fingerprint/data/"
	}

# load raw SCE and POE matrix, need to select appropriate dates
dataPath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/"

defineFile<-function(){
readline("enter filename: ")
}
dir(path = dataPath,
    pattern = "POE.matrix.20")
print("select POE file")
POEfile<-defineFile()

dir(path = dataPath,
pattern = ".20")
print("select SCE file")
SCEfile<-defineFile()
dir(path = dataPath,
pattern = "platform.frame.20")
print("select platform file")
platformfile<-defineFile()


load(paste(dataPath, POEfile, sep = ""))
load(paste(dataPath, SCEfile, sep = ""))
load(paste(dataPath, platformfile, sep = ""))

# N.B. the POE.matrix is slightly smaller as a few arrays were removed
dim(SCE.frame)
dim(POE.matrix)
dim(platform.frame)

# re-organize
SCE.frame<-t(SCE.frame)
colnames(SCE.frame)<-gsub("SCE_", "", colnames(SCE.frame))

platform.frame$GEO<-gsub("SCE_", "", platform.frame$GEO)
rownames(platform.frame)<-platform.frame$GEO

platform.frame<-platform.frame[rownames(platform.frame) %in% colnames(POE.matrix),]

SCE.frame<-SCE.frame[,colnames(SCE.frame) %in% colnames(POE.matrix)]

# check dimensions and orders
dim(SCE.frame)
dim(POE.matrix)
dim(platform.frame)
sum(colnames(SCE.frame) == colnames(POE.matrix))
sum(platform.frame$GEO == colnames(POE.matrix))

# now need to define thresholds for each platform
# N.B. this produces warnings and NaN values where pathway scores are NA
# That's okay for now
platforms<-unique(as.character(platform.frame$Platform))
# set threshold as 0.001
threshold<-0.001
platform.thresholds<-vector("list", length(platforms))
names(platform.thresholds)<-platforms
print("Check the that outputs are all close to defined threshold level")
for (i in 1:length(platform.thresholds)){
  # define submatrix for each platform
  POE.subMatrix<-POE.matrix[,platform.frame$Platform == names(platform.thresholds)[i]]
  SCE.subMatrix<-SCE.frame[,platform.frame$Platform == names(platform.thresholds)[i]]
  # setup matrix to contain threshold values for each pathway
  thresholdTable<-matrix(nrow = nrow(POE.subMatrix), ncol = 2)
  rownames(thresholdTable)<-rownames(POE.subMatrix)
  colnames(thresholdTable)<-c("high", "low")

# find closest reference to threshold
  for (j in 1:nrow(POE.subMatrix)){
    # obtain locations of points closest to chosen threshold
    highThresholdRef<-which(abs(threshold - POE.subMatrix[j,]) == min(abs(threshold - POE.subMatrix[j,]), na.rm = TRUE))
    lowThresholdRef<-which(abs((-threshold) - POE.subMatrix[j,]) == min(abs((-threshold) - POE.subMatrix[j,]), na.rm = TRUE))
    # use the mean value of these points as the threshold
    # print for error checking
    print(c(mean(POE.subMatrix[j,highThresholdRef]), mean(POE.subMatrix[j,highThresholdRef])))
    # insert into matrix
    thresholdTable[j,1]<-mean(SCE.subMatrix[j,highThresholdRef])
    thresholdTable[j,2]<-mean(SCE.subMatrix[j,lowThresholdRef])
  }
  platform.thresholds[[i]]<-thresholdTable
  }
  
save(platform.thresholds, file = paste(dataPath, "POE_derived_thresholds_0.001.",
    Sys.Date(), ".RData",
    sep = ""
    )
)

# also save into fingerprinting repository
try(system(paste("cp ",
             pathprintRepository,
             "platform.thresholds.RData ",
             dataPath,
             "platform.thresholds.RData.old",
             sep = "")))
try(system(paste("cp ",
             pathprintRepository,
             "GEO.fingerprint.matrix.RData ",
             dataPath,
             "GEO.fingerprint.matrix.RData.old",
             sep = "")))
save(platform.thresholds, file = paste(pathprintRepository, "platform.thresholds.RData", sep = ""))

GEO.fingerprint.matrix<-(POE.matrix > 0.001) - (POE.matrix < -0.001)
save(GEO.fingerprint.matrix, file = paste(pathprintRepository, "GEO.fingerprint.matrix.RData", sep = ""))

# Not necessary to do this last bit
#################
# save thresholded matrix and POE matrix into shared directory
# use thresholds 0.01 and 0.001
#POE.matrix.0.01<-(POE.matrix > 0.01) - (POE.matrix < -0.01)
#POE.matrix.0.001<-(POE.matrix > 0.001) - (POE.matrix < -0.001)
#POE.matrix.0.0001<-(POE.matrix > 0.0001) - (POE.matrix < -0.0001)

# how many zeros?
#(sum((POE.matrix.0.01 == 0), na.rm = TRUE))/(491 * 159586)
#(sum((POE.matrix.0.001 == 0), na.rm = TRUE))/(491 * 159586)
#(sum((POE.matrix.0.0001 == 0), na.rm = TRUE))/(491 * 159586)

#save(POE.matrix.0.01, file = "/data/shared/Fingerprint/POE.matrix.0.01.2011-04-06")
#save(POE.matrix.0.001, file = "/data/shared/Fingerprint/POE.matrix.0.001.2011-04-06")
#save(POE.matrix.0.0001, file = "/data/shared/Fingerprint/POE.matrix.0.0001.2011-04-06")

# END
