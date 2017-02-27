# Analysis of raw SCG matrix on the server
# source POE analysis package

load("/home/galtschu2/Documents/Projects/Fingerprinting/data/platform_SCG_frames.RData")


source("http://bioconductor.org/biocLite.R")
biocLite("metaArray", lib = "/home/galtschu2/R/x86_64-redhat-linux-gnu-library/2.11")

library(metaArray)


setwd("/home/galtschu2/Documents/Projects/Fingerprinting/")

fingerpath<-"/home/galtschu2/Documents/Projects/Fingerprinting/data/Fingerprints"
setwd(fingerpath)

files<-dir(fingerpath)[file.info(dir(fingerpath))$size > 200]
files<-files[grep("KEGG_Wiki_static_", files)]

# don't load langmoen files

files<-files[-(grep("LangmoenLangmoen", files))]



SCG<-vector("list", length(files))
platform<-vector("list", length(files))
pb <- txtProgressBar(min = 0, max = length((files)), style = 3)
for (i in 1:length(files)){
  load(files[i])
	setTxtProgressBar(pb, i)
	SCG[i]<-temp1[[1]]$SCG
  platform[i]<-temp1[[1]]$platform
  }

for (i in 1:5){
load(files[i])
setTxtProgressBar(pb, i)
SCG[i]<-temp1[[1]]$SCG
platform[i]<-temp1[[1]]$platform
}


if (!(exists("kegg_wiki_TR_static.Hs.gs"))){
load("/home/galtschu2/Documents/Databases/Gene sets/kegg_wiki_TR_static.Hs.gs", .GlobalEnv)
}


names(SCG)<-gsub("KEGG_Wiki_static_", "", files)
SCG.frame<-t(as.data.frame(SCG))
colnames(SCG.frame)<-names(kegg_wiki_TR_static.Hs.gs)

names(platform)<-gsub("KEGG_Wiki_static_", "", files)
platform.frame<-data.frame(GEO = names(platform), Platform = unlist(platform))
save(platform.frame, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/platform.frame.RData")
# top platforms are GPL570, GPL96, GPL1261, GPL81
table(platform.frame$Platform)

# plot histograms for top platforms, 
setwd("/home/galtschu2/Documents/Projects/Fingerprinting/data/")
pdf(file = "SCG.histgrams.pdf")
par(mfcol = c(2,2))
for (i in 1:ncol(SCG.frame)){
  hist(SCG.frame[grep("GPL570", platform.frame$Platform), i], nclass = 100,
        main = paste("GPL570", colnames(SCG.frame)[i]))
  hist(SCG.frame[grep("GPL96", platform.frame$Platform), i], nclass = 100,
        main = paste("GPL96", colnames(SCG.frame)[i]))
  hist(SCG.frame[grep("GPL1261", platform.frame$Platform), i], nclass = 100,
        main = paste("GPL1261", colnames(SCG.frame)[i]))
  hist(SCG.frame[grep("GPL81", platform.frame$Platform), i], nclass = 100,
        main = paste("GPL81", colnames(SCG.frame)[i]))
  }
dev.off()
  
# try calculating the POE (probability of expression) using the metaArray package
pb <- txtProgressBar(min = 0, max = ncol(SCG.frame), style = 3)

for (i in 1:ncol(SCG.frame)){
#for (i in 1:5){
  png(file = paste("SCG.POE.", i, ".png", sep = "")) # pdf 
  em.draw(t(SCG.frame[grep("GPL570", platform.frame$Platform), i]),
          cl = rep(0,length(grep("GPL570", platform.frame$Platform)))
          )
  setTxtProgressBar(pb, i)        
  dev.off()
  }

# ran this up to around 80 pathways then got bored - enough to get the general idea

# repeat, now calculating the POE for each pathway
SCG.POE<-matrix(nrow = nrow(SCG.frame), ncol = ncol(SCG.frame))
colnames(SCG.POE)<-colnames(SCG.frame)
rownames(SCG.POE)<-rownames(SCG.frame)
POE.Pi<-vector("numeric", nrow(SCG.frame))
names(POE.Pi)<-rownames(SCG.frame)

pb <- txtProgressBar(min = 0, max = ncol(SCG.frame), style = 3)
for (i in 1:ncol(SCG.frame)){
#for (i in 1:2){
  sample <- SCG.frame[grep("GPL570", platform.frame$Platform), i]
  fit<-fit.em(t(sample),
          cl = rep(0,length(sample))
          )
  SCG.POE[grep("GPL570", platform.frame$Platform), i]<-t(fit$expr)
  POE.Pi[grep("GPL570", platform.frame$Platform)][i]<-fit$Pi
  setTxtProgressBar(pb, i)        
  }


# save dataframe in case things go wrong from here on
save(SCG.POE, file = "SCG.POE.temp.RData")
save(POE.Pi, file = "POE.Pi.temp.RData")


# remaining arrays
remainingArrays<-as.character(unique(platform.frame$Platform))
remainingArrays<-remainingArrays[-grep("GPL570", remainingArrays)]

pb <- txtProgressBar(min = 0, max = ncol(SCG.frame), style = 3)
for (j in 1:length(remainingArrays)){
  array<-remainingArrays[j]
  print(paste("calculating POE for ", array, sep = ""))
  for (i in 1:ncol(SCG.frame)){
    sample <- SCG.frame[grep(array, platform.frame$Platform), i]
    fit<-fit.em(t(sample),
            cl = rep(0,length(sample))
            )
    SCG.POE[grep(array, platform.frame$Platform), i]<-t(fit$expr)
    POE.Pi[grep(array, platform.frame$Platform)][i]<-fit$Pi
    setTxtProgressBar(pb, i)        
    }
  }

save(SCG.POE, file = "SCG.POE.RData")
save(POE.Pi, file = "POE.Pi.RData")


# compare ecdf of full ditribution with the model
ecdf.3<-ecdf(SCG.frame[grep("GPL570", platform.frame$Platform), 3])

jpeg("test.plots.jpg")
par(mfcol = c(2,2))
hist(SCG.frame[grep("GPL570", platform.frame$Platform), 3], main = "Raw", nclass = 100, xlab = "Raw")
plot(SCG.frame[grep("GPL570", platform.frame$Platform), 3],
SCG.POE[grep("GPL570", platform.frame$Platform), 3], xlab = "Raw", ylab = "POE")
plot(ecdf(SCG.frame[grep("GPL570", platform.frame$Platform), 3]), main = "", xlab = "Raw", ylab = "ECDF")
plot(ecdf.3(SCG.frame[grep("GPL570", platform.frame$Platform), 3]),
SCG.POE[grep("GPL570", platform.frame$Platform), 3], xlab = "ECDF", ylab = "POE")

dev.off()

#####
# plot graphs

# plot POE thresholded data (p = 0.5)
pdf(file = "SCG.POE.histgrams.pdf")
par(mfcol = c(2,2))
for (i in 1:ncol(SCG.frame)){
#for (i in 1:10){
  platform.sub<-c("GPL570", "GPL96", "GPL1261", "GPL81")
  for (j in 1:4){
  data<-SCG.frame[grep(platform.sub[j], platform.frame$Platform), i]
  POE.data<-SCG.POE[grep(platform.sub[j], platform.frame$Platform), i]
  temp.hist<-hist(data,
                  nclass = 50,
                  main = paste(platform.sub[j], colnames(SCG.frame)[i])
                  )
             hist(data[POE.data>0.5],
                  breaks = temp.hist$breaks,
                  col = "red", add = TRUE
                  )
             hist(data[POE.data< c(-0.5)],
                  breaks = temp.hist$breaks,
                  col = "blue", add = TRUE
                  )                  
  }
  }
dev.off()

# plot quantile thresholded data (15%)

pdf(file = "SCG.quantile.histgrams.pdf")
par(mfcol = c(2,2))
for (i in 1:ncol(SCG.frame)){
#for (i in 1:10){
  platform.sub<-c("GPL570", "GPL96", "GPL1261", "GPL81")
  for (j in 1:4){
  data<-SCG.frame[grep(platform.sub[j], platform.frame$Platform), i]
  temp.hist<-hist(data,
                  nclass = 50,
                  main = paste(platform.sub[j], colnames(SCG.frame)[i])
                  )
             hist(data[data>quantile(data, 0.85)],
                  breaks = temp.hist$breaks,
                  col = "red", add = TRUE
                  )
             hist(data[data<quantile(data, 0.15)],
                  breaks = temp.hist$breaks,
                  col = "blue", add = TRUE
                  )                  
  }
  }
dev.off()

# only the top bit saved
save.image(file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/platform_SCG_frames.RData")



