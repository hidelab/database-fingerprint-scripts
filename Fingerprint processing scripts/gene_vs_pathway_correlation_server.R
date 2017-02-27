# plotting geneset gene-level vs pathway correlation
# Use pearson correlation or rank correlation..?
# Use thresholded or continuous matrix
# one platform or many platforms?
# how to plot - need upper-left diagonal

# testing
M<-matrix(sample(1:100,30), ncol = 3)
N<-matrix(sample(1:100,30), ncol = 3)
plot(cor(M), cor(N))
# this works fine
# now on server
# for now just compare correlation matrix for GPL570
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/platform_SCG_frames.RData")
load("/home/galtschu2/Documents/Projects/Ranking/data/GPL570.rank.matrix.RData")

# remove full fingerprint list from the workspace
rm(SCG)

GPL570.SCG<-t(SCG.frame[colnames(rank.matrix),])
save(GPL570.SCG, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/GPL570.SCG.RData")
q()
# quit and then re-open
load("/home/galtschu2/Documents/Projects/Ranking/data/GPL570.rank.matrix.RData")
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/GPL570.SCG.RData")
# trying to construct correlation matrices crashed first time around, attempt again but test 
# problem is memory limit - can't make a 45000x45000 matrix very easily!
# Try just the first 5000 arrays
# solve this by finding the correlation line by line and adding the points to the graph (png)
# i.e. 45000 iterations in a loop, with both rank and scg within loop
# Remeber to set-up graph for 0,1 then, plot in the x,y for row i etc
# from graph determine 'interesting' limits, then go back and extract these


rank.matrix.cor<-cor(rank.matrix[,1:5000])
GPL570.SCG.cor<-cor(GPL570.SCG[,1:5000])

png("/home/galtschu2/Documents/Projects/Fingerprinting/rank_V_SCG.png")
plot(rank.matrix.cor, GPL570.SCG.cor, xlab = "gene correlation", ylab = "pathway correlation")
dev.off()

# Looks like these are quite a few 
max(rank.matrix.cor[abs(GPL570.SCG.cor)<0.05])
max(GPL570.SCG.cor[abs(rank.matrix.cor)<0.05])

# Top zero correlation array is due to the misclassification of the species in GEO
which(rank.matrix.cor == max(rank.matrix.cor[abs(GPL570.SCG.cor)<0.05]) ,arr.ind = TRUE)

# this worked okay, now try to implement the larger dataset
# actually only need the upper diagonal elements of the correlation matrix
png("/home/galtschu2/Documents/Projects/Fingerprinting/rank_V_SCG_full.png")
# setup plot
par(pch = ".")
plot.new()
plot.window(xlim = c(-1,1), ylim = c(-1,1),
            xlab = "gene correlation", ylab = "pathway correlation")
axis(1)
axis(2)

# calculate correlation - just use 2 to start
pb<-txtProgressBar(min=0, max=2, style=3)
#for (i in 1:(length(ncol(rank.matrix))-1)){
for (i in 1:2){
  refcol<-(i+1):ncol(rank.matrix)
  cor.rank<-cor(rank.matrix[,i], rank.matrix[,refcol])
  cor.SCG<-cor(GPL570.SCG[,i], GPL570.SCG[,refcol])
  points(cor.rank[1,], cor.SCG[1,])
  setTxtProgressBar(pb, i)
}
dev.off()


#################
# Parallel does not work as cannot get script to write to same graph :(
# new idea, just output to a png and then combine using imagemagic
# testing running this in parallel using foreach
library(foreach)
library(multicore)
library(doMC)
registerDoMC(cores = 10)

foreach(i = 1:3) %dopar% sqrt(i)

ptime<-system.time({
  r <- foreach(icount(trials), .combine = cbind) %dopar% {
      ind <- sample(100, 100, replace = TRUE)
      result1 <- glm(x[ind, 2] ~ x[ind, 1], family = binomial(logit))
      coefficients(result1)
    }
  })[3]

stime<-system.time({
r <- foreach(icount(trials), .combine = cbind) %do% {
    ind <- sample(100, 100, replace = TRUE)
    result1 <- glm(x[ind, 2] ~ x[ind, 1], family = binomial(logit))
    coefficients(result1)
  }
})[3]

# runs in 13.453 vs 73.313s

# re-write script in parallel
png("/home/galtschu2/Documents/Projects/Fingerprinting/rank_V_SCG_full.png")
# setup plot
par(pch = ".")
plot.new()
plot.window(xlim = c(-1,1), ylim = c(-1,1),
            xlab = "gene correlation", ylab = "pathway correlation")
axis(1)
axis(2)

# calculate correlation - just use 2 to start
pb<-txtProgressBar(min=0, max=2, style=3)
#foreach (i = 1:(length(ncol(rank.matrix.cor))-1)){
foreach (i = 1:2) %dopar% {
  refcol<-(i+1):ncol(rank.matrix)
  cor.rank<-cor(rank.matrix[,i], rank.matrix[,refcol])
  cor.SCG<-cor(GPL570.SCG[,i], GPL570.SCG[,refcol])
  points(cor.rank[1,], cor.SCG[1,])
  setTxtProgressBar(pb, i)
}
dev.off()

png("/home/galtschu2/Documents/Projects/Fingerprinting/test.png")
# setup plot
par(pch = ".")
plot.new()
plot.window(xlim = c(-1,1), ylim = c(-1,1),
            xlab = "gene correlation", ylab = "pathway correlation")
axis(1)
axis(2)

pb<-txtProgressBar(min=0, max=2, style=3)
foreach (i = 1:10) %dopar% {
  eval(points(i/20, i/20), .GlobalEnv)
  print(dev.cur())
  setTxtProgressBar(pb, i)
}
dev.off()

plotPoints<-function(x){
  refcol<-(x+1):ncol(rank.matrix)
  cor.rank<-cor(rank.matrix[,x], rank.matrix[,refcol])
  cor.SCG<-cor(GPL570.SCG[,x], GPL570.SCG[,refcol])
  points(cor.rank[1,], cor.SCG[1,])
  setTxtProgressBar(pb, x)
}

testPar<-function(x){
    print(environment())
    eval(points(x/20, x/20), .GlobalEnv)
    eval(print(environment()), testenv)
    print(dev.cur())
    setTxtProgressBar(pb, x)
  }

testenv <- new.env()
png("/home/galtschu2/Documents/Projects/Fingerprinting/test.png")
# setup plot

plot.new()
plot.window(xlim = c(-1,1), ylim = c(-1,1),
            xlab = "gene correlation", ylab = "pathway correlation")
axis(1)
axis(2)
par(pch = ".", new = TRUE)
pb<-txtProgressBar(min=0, max=2, style=3)
mclapply(1:10, testPar)
dev.off()


########## This works okay now

pb<-txtProgressBar(min=0, max=2, style=3)
foreach (i = 1:20) %dopar% {
  png(paste("/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/rankVSCG", i, "png", sep = "."))
  par(bg = "transparent")
  plot(i/20, i/20, xlim = c(-1,1), ylim = c(-1,1))
  dev.off()
  setTxtProgressBar(pb, i)
}

pb<-txtProgressBar(min=0, max=10, style=3)
#foreach (i == 1:(length(ncol(rank.matrix))-1)){
foreach (i = 1:10) %dopar% {
  refcol<-(i+1):ncol(rank.matrix)
  cor.rank<-cor(rank.matrix[,i], rank.matrix[,refcol])
  cor.SCG<-cor(GPL570.SCG[,i], GPL570.SCG[,refcol])
  png(paste("/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/rankVSCG", i, "png", sep = "."))
  par(bg = "transparent")
  plot(cor.rank[1,], cor.SCG[1,])
  dev.off()
  setTxtProgressBar(pb, i)
}

# this is very memory intensive!
# try just using mclapply
multiPlots<-function(x){
  refcol<-(x+1):ncol(rank.matrix)
  cor.rank<-cor(rank.matrix[,x], rank.matrix[,refcol])
  cor.SCG<-cor(GPL570.SCG[,x], GPL570.SCG[,refcol])
  png(paste("/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/rankVSCG", x, "png", sep = "."))
  par(bg = "transparent")
  plot(cor.rank[1,], cor.SCG[1,])
  dev.off()
  setTxtProgressBar(pb, x)
}

mclapply(1:5, multiPlots)

# this doesn't seem to help
# perhaps there is another way around this?

# try to break matrix down into smaller pieces
# Try split by 10
n<-10
splits<-lapply(1:10, function(x){(1+(ncol(rank.matrix) %/% n)*(x-1)):((ncol(rank.matrix) %/% n)*x)})
splits[[n]]<-c(splits[[10]], tail(1:ncol(rank.matrix), ncol(rank.matrix) %% n))





rank.matrix.1<-rank.matrix[,splits[[1]]]
rank.matrix.2<-rank.matrix[,splits[[2]]]

GPL570.SCG.1<-GPL570.SCG[,splits[[1]]]
GPL570.SCG.2<-GPL570.SCG[,splits[[2]]]

pb<-txtProgressBar(min=0, max=10, style=3)
#foreach (i == 1:(length(ncol(rank.matrix))-1)){
system.time(foreach (i = 1:10) %dopar% {
  refcol<-(i+1):ncol(rank.matrix.1)
  cor.rank<-cor(rank.matrix.1[,i], rank.matrix.1[,refcol])
  cor.SCG<-cor(GPL570.SCG.1[,i], GPL570.SCG.1[,refcol])
  png(paste("/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/rankVSCG", i, "png", sep = "."))
  par(bg = "transparent", pch = ".")
  plot(cor.rank[1,], cor.SCG[1,])
  dev.off()
  setTxtProgressBar(pb, i)
})

multiPlots.1<-function(x){
  refcol<-(x+1):ncol(rank.matrix.1)
  cor.rank<-cor(rank.matrix.1[,x], rank.matrix.1[,refcol])
  cor.SCG<-cor(GPL570.SCG.1[,x], GPL570.SCG.1[,refcol])
  png(paste("/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/rankVSCG", x, i ,j ,"png", sep = "."))
  par(bg = "transparent", pch = ".")
  plot.new()
  plot.window(xlim = c(-1,1), ylim = c(-1,1),
              xlab = "gene correlation", ylab = "pathway correlation")
  axis(1)
  axis(2)
  points(cor.rank[1,], cor.SCG[1,])
  dev.off()
  setTxtProgressBar(pb, x)
}


system.time(mclapply(1:ncol(rank.matrix.1), multiPlots.1, mc.cores = 10))

# actually seems very similar times, need to specify number of cores with mclapply
# now need to set up a loop
i<-1
j<-1
pb<-txtProgressBar(min=0, max=ncol(rank.matrix.1), style=3)
mclapply(1:ncol(rank.matrix.1), multiPlots.1, mc.cores = 10)

# this seems to run okay
# now try to set up loop of loops

multiPlots.2<-function(x){
  refcol<-(x+1):ncol(rank.matrix.1)
  cor.rank<-cor(rank.matrix.1[,x], rank.matrix.1[,refcol])
  cor.SCG<-cor(GPL570.SCG.1[,x], GPL570.SCG.1[,refcol])
  png(paste("/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/rankVSCG", x, i ,j ,"png", sep = "."))
  par(bg = "transparent", pch = ".")
  plot.new()
  plot.window(xlim = c(-1,1), ylim = c(-1,1),
              xlab = "gene correlation", ylab = "pathway correlation")
  axis(1)
  axis(2)
  points(cor.rank[1,], cor.SCG[1,])
  dev.off()
  setTxtProgressBar(pb, x)
}

# additional adjustment for memory
# can make the objects all integers - max 10 decimal places
# determine the limits

object.size(rank.matrix)
object.size(GPL570.SCG)

max(rank.matrix)
min(rank.matrix)
max(GPL570.SCG)
min(GPL570.SCG)

# multiply by a suitable number and then turn into an integer
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/GPL570.SCG.RData")
load("/home/galtschu2/Documents/Projects/Ranking/data/GPL570.rank.matrix.RData")

should be ncol - 1
rank.matrix.int<-apply(1000*rank.matrix, c(1,2), as.integer)
GPL570.SCG.int <-apply(1000*GPL570.SCG, c(1,2), as.integer)
rank.matrix.int<-matrix(as.integer(1000*rank.matrix), ncol = ncol(rank.matrix))
GPL570.SCG.int<-matrix(as.integer(1000*GPL570.SCG), ncol = ncol(GPL570.SCG))
object.size(rank.matrix.int)
object.size(GPL570.SCG.int)

save(GPL570.SCG.int, file = "/home/galtschu2/Documents/Projects/Fingerprinting/data/GPL570.SCG.int.RData")
save(rank.matrix.int, file = "/home/galtschu2/Documents/Projects/Ranking/data/rank.matrix.int.RData")

# start again with integer matrices
load("/home/galtschu2/Documents/Projects/Fingerprinting/data/GPL570.SCG.int.RData")
load("/home/galtschu2/Documents/Projects/Ranking/data/rank.matrix.int.RData")

path<-"/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/"
setwd(path)

library(foreach)
library(multicore)
library(doMC)
registerDoMC(cores = 10)

multiPlots.diag<-function(x){
  refcol<-(x+1):ncol(rank.matrix.j)
  cor.rank<-cor(rank.matrix.i[,x], rank.matrix.j[,refcol])
  cor.SCG<-cor(GPL570.SCG.i[,x], GPL570.SCG.j[,refcol])
  png(paste("/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/rankVSCG", x, i ,j ,"png", sep = "."))
  par(bg = "transparent", pch = ".")
  plot.new()
  plot.window(xlim = c(-1,1), ylim = c(-1,1),
              xlab = "gene correlation", ylab = "pathway correlation")
  axis(1)
  axis(2)
  points(cor.rank[1,], cor.SCG[1,], col = "#00000050")
  dev.off()
  setTxtProgressBar(pb, x)
}

# corrected to allow for all combinations of i,j
# means double counting but that's too bad
multiPlots<-function(x){
  cor.rank<-cor(rank.matrix.i[,x], rank.matrix.j)
  cor.SCG<-cor(GPL570.SCG.i[,x], GPL570.SCG.j)
  png(paste("/home/galtschu2/Documents/Projects/Fingerprinting/correlationOutput/rankVSCG", x, i ,j ,"png", sep = "."))
  par(bg = "transparent", pch = ".")
  plot.new()
  plot.window(xlim = c(-1,1), ylim = c(-1,1),
              xlab = "gene correlation", ylab = "pathway correlation")
  axis(1)
  axis(2)
  points(cor.rank[1,], cor.SCG[1,], col = "#00000050")
  dev.off()
  setTxtProgressBar(pb, x)
}


n<-10
splits<-lapply(1:10, function(x){(1+(ncol(rank.matrix.int) %/% n)*(x-1)):((ncol(rank.matrix.int) %/% n)*x)})
splits[[n]]<-c(splits[[10]], tail(1:ncol(rank.matrix.int), ncol(rank.matrix.int) %% n))





pb.i<-txtProgressBar(min=0, max=n, style=2)
pb.j<-txtProgressBar(min=0, max=n, style=1)

for(i in 1){
#for(i in 1:n){
  for (j in 1){
  #for (j in 1:n){
    rank.matrix.i<-rank.matrix.int[,splits[[i]]]
    GPL570.SCG.i<-GPL570.SCG.int[,splits[[i]]]     
    rank.matrix.j<-rank.matrix.int[,splits[[j]]]
    GPL570.SCG.j<-GPL570.SCG.int[,splits[[j]]]
    pb<-txtProgressBar(min=0, max=(ncol(rank.matrix.i)), style=3)
    ##if (i == j){
      ##mclapply(1:(ncol(rank.matrix)-4350), multiPlots.diag, mc.cores = 4)
      ##mclapply(1:(ncol(rank.matrix)), multiPlots.diag, mc.cores = 4)
    ##  }
    ##else if(!(i == j)){
      #mclapply(1:(ncol(rank.matrix)-4350), multiPlots, mc.cores = 4)
      mclapply(1:(ncol(rank.matrix.i)), multiPlots, mc.cores = 4)
    ##  }
    #system(paste("convert -page rank* -flatten ./compiled/", i, ".", j, ".", "convert.png", sep = ""))
    #system("rm rank*")
    print(j)
    setTxtProgressBar(pb.j, j)
      }
    print(i)
    setTxtProgressBar(pb.i, i)
}

# bug hunting - for some reason zero pathway correlation arrays did not come up
# e.g. GSM132948

match("GSM132948", colnames(GPL570.SCG.int))
# no colnames!

# conversion to integer made all of the very small SCG values into zero
# correlation with zero vector gives NA so these graphs are not plotted
# hence the 'line' across at zero pathway correation does not come up
# this is a good as differentiates between zero pathawy correlation due to
# actual correlation and due to a zero vector

# 72.1.1 gives interesting pattern, mainly negative pathway correlation,
# zero gene correlation, with pathway correlation as high as 0.5 for some arrays

multiDiffCorrelation<-function(x){
  cor.diff<-cor(GPL570.SCG.i[,x], GPL570.SCG.j) - cor(rank.matrix.i[,x], rank.matrix.j)
  if(sum(cor.diff > 0.5, na.rm = TRUE)>0){print(paste(x, i, j, which(cor.diff > 0.5), sep = "."))}
  setTxtProgressBar(pb, x)
}

for(i in 1){
#for(i in 1:n){
  for (j in 1){
  #for (j in 1:n){
    rank.matrix.i<-rank.matrix.int[,splits[[i]]]
    GPL570.SCG.i<-GPL570.SCG.int[,splits[[i]]]     
    rank.matrix.j<-rank.matrix.int[,splits[[j]]]
    GPL570.SCG.j<-GPL570.SCG.int[,splits[[j]]]
    pb<-txtProgressBar(min=0, max=(ncol(rank.matrix.i)), style=3)
    ##if (i == j){
      ##mclapply(1:(ncol(rank.matrix)-4350), multiPlots.diag, mc.cores = 4)
      ##mclapply(1:(ncol(rank.matrix)), multiPlots.diag, mc.cores = 4)
    ##  }
    ##else if(!(i == j)){
      ##mclapply(1:(ncol(rank.matrix)-4350), multiPlots, mc.cores = 4)
      #mclapply(1:(ncol(rank.matrix.i)), multiDiffCorrelation, mc.cores = 2)
      lapply(1:(ncol(rank.matrix.i)), multiDiffCorrelation)
  ##  }
    #system(paste("convert -page rank* -flatten ./compiled/", i, ".", j, ".", "convert.png", sep = ""))
    #system("rm rank*")
    print(j)
    setTxtProgressBar(pb.j, j)
      }
    print(i)
    setTxtProgressBar(pb.i, i)
}


# examples of highly differential correlation
# 558.1.1.598"  "558.1.1.2403" "558.1.1.2407" "558.1.1.2408
# 561.1.1.598"  "561.1.1.2403" "561.1.1.2407" "561.1.1.2408
# 563.1.1.2403


cor(GPL570.SCG.int[,558], GPL570.SCG.int[,2403])
0.7696118
cor(rank.matrix.int[,558], rank.matrix.int[,2403])
0.2528505

cor(GPL570.SCG.int[,558], GPL570.SCG.int[,598])
0.7696118
cor(rank.matrix.int[,558], rank.matrix.int[,598])
0.1673820

cor(GPL570.SCG.int[,558], GPL570.SCG.int[,2407])
0.7391786
cor(rank.matrix.int[,558], rank.matrix.int[,2407])
0.2386444

colnames(GPL570.SCG)[2407]