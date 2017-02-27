GSE24759.dat <- phyDat(t(GSE24759.consensus), type = "USER", levels = c(-1,0,1))
GSE24759.parsimony<-pratchet(
  GSE24759.dat, start = GSE24759.NJ.tree, k = 50,
  method = "sankoff", cost = CM, trace = 0, np = 1, all = TRUE)
parsimony(GSE24759.parsimony, GSE24759.dat)


GSE7658.inform<-rownames(GSE7658.consensus)[apply(GSE7658.consensus, 1, sd) > 0]
GSE7658.inform1<-rownames(GSE7658.data)[apply(GSE7658.data, 1, sd) > 0]

GSE7658.top.tree<-list()
GSE7658.min<-list()
for (i in 1:6){
name<-colnames(GSE7658.data)[i]
#dat.add<-phyDat(t(cbind(GSE24759.consensus, GSE7658.data[,i,drop = FALSE])[GSE24759.inform,]), type = "USER", levels = c(-1,0,1))
dat.add<-phyDat(t(cbind(GSE24759.consensus, GSE7658.data[,i,drop = FALSE])), type = "USER", levels = c(-1,0,1))
invisible(GSE24759.parsimony.add<-add.everywhere(GSE24759.parsimony, name))
p<-parsimony(GSE24759.parsimony.add, dat.add, method = "sankoff", cost = CM)
min[[i]]<-which(p == min(p))
  for (j in 1:length(min[[i]])){
  GSE7658.top.tree<-append(GSE7658.top.tree, list(GSE24759.parsimony.add[[min[[i]][j]]]))
    }
  }

GSE7658.top.tree<-lapply(GSE7658.top.tree, root, 21, resolve.root = TRUE)
class(GSE7658.top.tree)<-"multiPhylo"
plot(GSE7658.top.tree, tip.color = c(rep("black", 38), "red"))


colnames(GSE6506.data)<-paste(colnames(GSE6506.data), GSE6506.meta$cellType)

GSE6506.top.tree<-list()
GSE6506.min<-list()
for (i in 1:10){
name<-colnames(GSE6506.consensus)[i]
dat.add<-phyDat(t(cbind(GSE24759.consensus, GSE6506.consensus[,i,drop = FALSE])[intersect(GSE6506.inform, GSE24759.inform),]), type = "USER", levels = c(-1,0,1))
#dat.add<-phyDat(t(cbind(GSE24759.consensus, GSE6506.data[,i,drop = FALSE])), type = "USER", levels = c(-1,0,1))
invisible(GSE24759.parsimony.add<-add.everywhere(GSE24759.parsimony, name))
p<-parsimony(GSE24759.parsimony.add, dat.add, method = "sankoff", cost = CM)
min[[i]]<-which(p == min(p))
  for (j in 1:length(min[[i]])){
  GSE6506.top.tree<-append(GSE6506.top.tree, list(GSE24759.parsimony.add[[min[[i]][j]]]))
    }
  }

GSE6506.top.tree<-lapply(GSE6506.top.tree, root, 21, resolve.root = TRUE)
class(GSE6506.top.tree)<-"multiPhylo"
plot(GSE6506.top.tree, tip.color = c(rep("black", 38), "red"))


setwd("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/Blood")
bloodSet<-read.delim("data/bloodTypesHumanMouse.txt", stringsAsFactors = FALSE)

bloodSet$name<-paste(bloodSet$type, bloodSet$species)
bloodSet.consensus<-sapply(unique(bloodSet$name), function(x){
    consensusFingerprint(GEO.fingerprint.matrix[,
        bloodSet$gsm[bloodSet$name == x]],
    threshold = threshold)
    })

bloodSet.human<-bloodSet.consensus[,grep("human", colnames(bloodSet.consensus))]
bloodSet.human.dat <- phyDat(t(bloodSet.human), type = "USER", levels = c(-1,0,1))
bloodSet.human.dist <- dist.hamming(bloodSet.human.dat)
# construct trees
bloodSet.human.NJ.tree <- NJ(bloodSet.human.dist)
bloodSet.human.parsimony <-pratchet(
  bloodSet.human.dat, start = bloodSet.human.NJ.tree, k = 50,
  method = "sankoff", cost = CM, trace = 0, np = 1, all = TRUE)


bloodSet.human.inform<-rownames(bloodSet.human)[apply(bloodSet.human, 1, sd) > 0]

# comparing mouse samples to mouse tree
GSE6506.top.tree<-list()
GSE6506.min<-list()
for (i in 1:10){
name<-colnames(GSE6506.consensus)[i]
#dat.add<-phyDat(t(cbind(bloodSet.human, GSE6506.consensus[,i,drop = FALSE])[intersect(GSE6506.inform, bloodSet.human.inform),]), type = "USER", levels = c(-1,0,1))
dat.add<-phyDat(t(cbind(bloodSet.human, GSE6506.consensus[,i,drop = FALSE])), type = "USER", levels = c(-1,0,1))
invisible(bloodSet.human.parsimony.add<-add.everywhere(bloodSet.human.parsimony, name))
p<-parsimony(bloodSet.human.parsimony.add, dat.add, method = "sankoff", cost = CM)
min[[i]]<-which(p == min(p))
  for (j in 1:length(min[[i]])){
  GSE6506.top.tree<-append(GSE6506.top.tree, list(bloodSet.human.parsimony.add[[min[[i]][j]]]))
    }
  }

GSE6506.top.tree<-lapply(GSE6506.top.tree, root, 4, resolve.root = TRUE)
class(GSE6506.top.tree)<-"multiPhylo"
plot(GSE6506.top.tree, tip.color = c(rep("black", 10), "red"))


GSE7658.data.inform<-rownames(GSE7658.data)[apply(GSE7658.data, 1, sd) > 0]
# comparing zfish samples to mouse tree
GSE7658.top.tree<-list()
GSE7658.min<-list()
for (i in 1:6){
name<-colnames(GSE7658.data)[i]
dat.add<-phyDat(t(cbind(bloodSet.human, GSE7658.data[,i,drop = FALSE])[intersect(bloodSet.human.inform, GSE7658.data.inform),]), type = "USER", levels = c(-1,0,1))
#dat.add<-phyDat(t(cbind(bloodSet.human, GSE7658.data[,i,drop = FALSE])), type = "USER", levels = c(-1,0,1))
invisible(bloodSet.human.parsimony.add<-add.everywhere(bloodSet.human.parsimony, name))
p<-parsimony(bloodSet.human.parsimony.add, dat.add, method = "sankoff", cost = CM)
min[[i]]<-which(p == min(p))
  for (j in 1:length(min[[i]])){
  GSE7658.top.tree<-append(GSE7658.top.tree, list(bloodSet.human.parsimony.add[[min[[i]][j]]]))
    }
  }

GSE7658.top.tree<-lapply(GSE7658.top.tree, root, 4, resolve.root = TRUE)
class(GSE7658.top.tree)<-"multiPhylo"
plot(GSE7658.top.tree, tip.color = c(rep("black", 10), "red"))


The above both provide very interesting examples of how the human and mouse trees can be combined. The problem will be finding a way to viz them














GSE6506.consensus[1:5, 1:5]




tree<-GSE24759.parsimony
tip.name<-"temp"


for(i in 1:3){
       trees[[i]]<-bind.tree(tree,new.tip,where=tree$edge[i,2],
  position=0.5)
       trees[[i]]$edge.length<-NULL
    }

test<-bind.tree(tree,new.tip,where=i, position=0.5); i<-i+1
bind.tree(tree,new.tip,where=30, position=0.5)
trees<-vector("list", 74); class(trees)<-"multiPhylo"

for (i in 1:74){
try(trees[[i]]<-bind.tree(tree,new.tip,where=i, position=0.5))
}

add.everywhere<-function(tree,tip.name){
  if(!require(ape)) stop("function needs 'ape' package.")
	if(class(tree)!="phylo") stop("tree should be an object of class 'phylo.'")
	tree<-reorder(tree, "cladewise") # convert to cladewise
  tree<-unroot(tree) # unroot tree
	tree$edge.length<-rep(1,nrow(tree$edge)) # set all edge lengths to 1.0
  # create new tip
	new.tip<-list(edge=matrix(c(2L,1L),1,2),tip.label=tip.name,edge.length=1,Nnode=1L)
	class(new.tip)<-"phylo"
  # add the new tip to all edges of the tree
	trees<-list(); class(trees)<-"multiPhylo"
	for(i in 1:nrow(tree$edge)){
		#try(trees[[i]]<-bind.tree(read.tree(text = write.tree(tree, file = "")),new.tip,where=tree$edge[i,2],position=0.5))
    try(trees[[i]]<-bind.tree(tree,new.tip,where=tree$edge[i,2],position=0.5))
		try(trees[[i]]$edge.length<-NULL)
	}
	trees<-trees[sapply(trees, class) == "phylo"]
  return(trees)
}

GSE24759$type[grep("CD4+ Effector Memory", GSE24759$type, fixed = TRUE)]<-"CD4+_activated"
GSE24759$type[setdiff(grep("CD8+ Effector Memory", GSE24759$type, fixed = TRUE),
                      grep("CD8+ Effector Memory RA", GSE24759$type, fixed = TRUE))
              ]<-"CD8+_activated"
GSE24759$type[grep("Naive CD4+ T-cell", GSE24759$type, fixed = TRUE)]<-"CD4+_naive"
GSE24759$type[grep("Naive CD8+ T-cell", GSE24759$type, fixed = TRUE)]<-"CD8+_naive"
GSE24759$type[grep("Mature B-cells", GSE24759$type, fixed = TRUE)]<-"B-Cell"
GSE24759$type[grep("Mature NK cell", GSE24759$type, fixed = TRUE)]<-"NK"
GSE24759$type[grep("Hematopoietic stem cell_CD133+", GSE24759$type, fixed = TRUE)]<-"LT_HSC"
GSE24759$type[setdiff(grep("Monocyte", GSE24759$type, fixed = TRUE),
                      grep("Colony Forming Unit-Monocyte", GSE24759$type, fixed = TRUE))
              ]<-"Monocyte"
GSE24759$type[grep("Granulocyte (Neutrophil)", GSE24759$type, fixed = TRUE)]<-"Granulocyte"
GSE24759$type[grep("Erythroid_CD34- CD71+", GSE24759$type, fixed = TRUE)]<-"Nucleated Erythrocytes"

GSE6506$type<-c(rep("LT_HSC", 2),
                rep("NK",2),
                rep("CD4+_naive", 2),
                rep("CD8+_naive", 2),
                rep("CD4+_activated", 2),
                rep("CD8+_activated", 2),
                rep("B-Cell", 2),
                rep("Monocyte", 2),
                rep("Granulocyte", 2),
                rep("Nucleated Erythrocytes", 2)
                )





# stuff removed from Sweave doc
GSE7658.consensus<-sapply(GSE7658.cellType_main, function(x){
    consensusFingerprint(GEO.fingerprint.matrix[,
        GSE7658.meta$GSM[GSE7658.meta$cellType_main == x]],
    threshold = threshold)
    })
GSE7658.consensus[1:5,]

GSE24759.inform<-rownames(GSE24759.consensus)[apply(GSE24759.consensus, 1, sd) > 0]
GSE7658.inform<-rownames(GSE7658.consensus)[apply(GSE7658.consensus, 1, sd) > 0]






  
for (i in 1:length(GSE7658.branches.inform)){
  plot(GSE7658.branches.inform[[i]], tip.color = GSE7658.branches.inform[[i]]$tip.color)
}
  
  
  
GSE7658.top.tree<-list()
GSE7658.min<-list()
for (i in ncol(GSE7658.data)){
name<-colnames(GSE7658.data)[i]
dat.add<-phyDat(t(cbind(bloodSet.human, GSE7658.data[,i,drop = FALSE])[GSE7658.data.inform,]), type = "USER", levels = c(-1,0,1))
#dat.add<-phyDat(t(cbind(bloodSet.human, GSE7658.data[,i,drop = FALSE])), type = "USER", levels = c(-1,0,1))
invisible(bloodSet.human.parsimony.add<-add.everywhere(bloodSet.human.parsimony, name))
p<-parsimony(bloodSet.human.parsimony.add, dat.add, method = "sankoff", cost = CM)
min[[i]]<-which(p == min(p))
  for (j in 1:length(min[[i]])){
  GSE7658.top.tree<-append(GSE7658.top.tree, list(bloodSet.human.parsimony.add[[min[[i]][j]]]))
    }
  }

GSE7658.top.tree<-lapply(GSE7658.top.tree, root, 4, resolve.root = TRUE)
class(GSE7658.top.tree)<-"multiPhylo"
plot(GSE7658.top.tree, tip.color = c(rep("black", 10), "red"))




Alternative approach is to add a sample one-by-one into the tree - this seems to give similar results
NEXT STEP IS TO ATTEMPT TO INTEGRATE THE MOUSE DATA IN THE SAME WAY - A ONE-BY-ONE APPROACH MAY BE BENEFICIAL TO PREVENT DISRUPTION OF THE HUMAN TREE
<<one_by_one>>=
for (i in ncol(GSE7658.consensus)){
  dat<-phyDat(t(cbind(GSE7658.consensus[,i,drop = FALSE], GSE24759.consensus)[intersect(GSE24759.inform,GSE7658.inform),]), type = "USER", levels = c(-1,0,1))
  dist <- dist.hamming(dat)
  NJ.tree <- NJ(dist)
  parsimony <- pratchet(
   dat, start = NJ.tree, k = 50,
    method = "sankoff", cost = CM, trace = 0, np = 1, all = TRUE)
}
@
<<GSE6506>>=
GSE6506.meta<-GEO.metadata.matrix[
    GEO.metadata.matrix$GSE %in% "GSE6506",]
GSE6506.meta$cellType<-sapply(GSE6506.meta$Source,
            function(x){unlist(strsplit(x, split = " isolated"))[[1]]})
GSE6506.meta$cellType<-gsub(" activated with an LPS treatment and", "", GSE6506.meta$cellType)
GSE6506.cellTypes<-levels(as.factor(GSE6506.meta$cellType))
@
The fingerprints can be extracted from the fingerprint matrix and a consensus fingerprint constructed for each of the cell types.
<<GSE6506.data>>=
GSE6506.data<-GEO.fingerprint.matrix[,GSE6506.meta$GSM]
GSE6506.consensus<-sapply(GSE6506.cellTypes, function(x){
    consensusFingerprint(GEO.fingerprint.matrix[,
        GSE6506.meta$GSM[GSE6506.meta$cellType == x]],
    threshold = threshold)
    })
GSE6506.consensus[1:5, 1:5]
@ 
<<GSE6506_fit>>=
GSE6506.inform<-rownames(GSE6506.consensus)[apply(GSE6506.consensus, 1, sd) > 0]
for (i in ncol(GSE6506.consensus)){
  dat<-phyDat(t(cbind(GSE6506.consensus[,i,drop = FALSE], GSE24759.consensus)[intersect(GSE24759.inform,GSE6506.inform),]), type = "USER", levels = c(-1,0,1))
  dist <- dist.hamming(dat)
  NJ.tree <- NJ(dist)
  parsimony <- pratchet(
   dat, start = NJ.tree, k = 50,
    method = "sankoff", cost = CM, trace = 0, np = 1, all = TRUE)
}

Not quite finished with this bit above
  
  

<<combined>>=
combined.data<-cbind(GSE7658.consensus, GSE24759.consensus)[intersect(GSE24759.inform,GSE7658.inform),]
combined.dat <- phyDat(t(combined.data), type = "USER", levels = c(-1,0,1))
combined.dist <- dist.hamming(combined.dat)
# construct trees
combined.NJ.tree <- NJ(combined.dist)
combined.NJ.tree<-root(combined.NJ.tree, 23, resolve.root = TRUE)
plot(combined.NJ.tree)

combined.parsimony <- pratchet(
  combined.dat, start = combined.NJ.tree, k = 50,
  method = "sankoff", cost = CM, trace = 0, np = 1, all = TRUE)

combined.parsimony<-root(combined.parsimony, 29, resolve.root = TRUE)

@




Can these blood types be matched to the GEO corpus? First try with the HSCs.
<<matchTypes>>=
# invisible(GSE24759.types.distance<-apply(GSE24759.consensus, 2,
#         consensusDistance, GEO.fingerprint.matrix))
# 
# GSE24759.types.distance.meta<-lapply(
#   names(GSE24759.types.distance), function(x){
#     cbind(GSE24759.types.distance[[x]],
#             GEO.metadata.matrix[
#               match(rownames(GSE24759.types.distance[[x]]),
#               GEO.metadata.matrix$GSM)
#                 ,])})
# names(GSE24759.types.distance.meta)<-names(GSE24759.types.distance)
@