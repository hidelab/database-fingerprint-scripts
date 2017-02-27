

<<brainQuantification>>=
sum(brain.meta$pvalue < 0.02)
sum(brain.meta$distance < 0.2)
top.brain<-brain.meta[brain.meta$pvalue < 0.02,]
top.brain<-brain.meta[1:7500,]


brain.valid<-(rowSums(
  sapply(brain.terms, function(x){
    grepl(x, apply(
            top.brain[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
brain.valid[top.brain$GSE %in% brain.GSE]<-1
non.brain<-top.brain[brain.valid == 0,]
View(non.brain)
sum(brain.valid)/length(brain.valid)
table(top.brain[brain.valid == 1,"Species"])
table(top.brain[brain.valid == 1,"GPL"])
gpl.name[match(names(table(top.brain[brain.valid == 1,"GPL"])), gpl.name$gpl),]
@
<<liverQuantification>>=
sum(liver.meta$pvalue < 0.02)
top.liver<-liver.meta[liver.meta$pvalue < 0.02,]
top.liver<-liver.meta[1:7500,]


# N.B. could use GSE12189 as a positive/negative control for liver

liver.valid<-(rowSums(
  sapply(liver.terms, function(x){
    grepl(x, apply(
            top.liver[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
liver.valid[top.liver$GSE %in% liver.GSE]<-1
non.liver<-top.liver[liver.valid == 0,]
View(non.liver)
sum(liver.valid)/length(liver.valid)
table(top.liver[liver.valid == 1,"Species"])
table(top.liver[liver.valid == 1,"GPL"])
gpl.name[match(names(table(top.liver[liver.valid == 1,"GPL"])), gpl.name$gpl),]
@
<<lungQuantification>>=
lung.meta<-cbind(tissueDistance$lung,
              GEO.metadata.matrix[
                match(rownames(tissueDistance$lung),
                GEO.metadata.matrix$GSM)
              ,])
sum(lung.meta$pvalue < 0.02)
sum(lung.meta$distance < 0.2)

top.lung<-lung.meta[lung.meta$pvalue < 0.02,]
top.lung<-lung.meta[lung.meta$distance < 0.15,]
top.lung<-lung.meta[1:7500,]

non.lung.terms<-c("breast", "bone marrow", "pancreas", "Ovarian", "colon")
# N.B. could use GSE12189 as a positive/negative control for lung
lung.valid<-(rowSums(
  sapply(lung.terms, function(x){
    grepl(x, apply(
            top.lung[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
lung.non.valid<-(rowSums(
  sapply(non.lung.terms, function(x){
    grepl(x, apply(
            top.lung[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)



lung.valid[top.lung$GSE %in% lung.GSE]<-1
non.lung<-top.lung[lung.valid == 0,]
View(top.lung[(lung.valid == 0 & lung.non.valid == 0),])
# checked top 500 out of 3679
View(non.lung)

sum(lung.valid)/length(lung.valid)
table(top.lung[lung.valid == 1,"Species"])
table(top.lung[lung.valid == 1,"GPL"])
gpl.name[match(names(table(top.lung[lung.valid == 1,"GPL"])), gpl.name$gpl),]
@
<<spleenQuantification>>=
spleen.meta<-cbind(tissueDistance$spleen,
              GEO.metadata.matrix[
                match(rownames(tissueDistance$spleen),
                GEO.metadata.matrix$GSM)
              ,])
sum(spleen.meta$pvalue < 0.02)
sum(spleen.meta$distance < 0.2)
top.spleen<-spleen.meta[1:7500,]
non.spleen.terms<-c("breast", "pancreas", "Ovarian", "colon", "lung", "Pulmonary")
spleen.blood.terms<-c(spleen.terms, "blood", "leukocytes", "neutrophils", "monocytes", "Dendritic",
                      "neutrophil", "PBMC", "macrophages", "reticulocyte", "PMN_Wegener's",
                      "monopcytes", "Buffycoat", "mononuclear", "leukemia", "Hematopoietic")
spleen.blood.GSE<-c(spleen.GSE, "GSE19743", "GSE7116", "GSE11907", "GSE2350", "GSE9960", "GSE22630")
# N.B. could use GSE12189 as a positive/negative control for spleen
spleen.valid<-(rowSums(
  sapply(spleen.terms, function(x){
    grepl(x, apply(
            top.spleen[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
spleen.non.valid<-(rowSums(
  sapply(non.spleen.terms, function(x){
    grepl(x, apply(
            top.spleen[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
spleen.blood.valid<-(rowSums(
  sapply(spleen.blood.terms, function(x){
    grepl(x, apply(
            top.spleen[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)            
             
             
spleen.valid[top.spleen$GSE %in% spleen.GSE]<-1
spleen.blood.valid[top.spleen$GSE %in% spleen.blood.GSE]<-1
non.spleen<-top.spleen[spleen.valid == 0,]
non.blood.spleen<-top.spleen[spleen.blood.valid == 0,]
View(non.spleen)
View(top.spleen[(spleen.valid == 0 & spleen.non.valid == 0),])
View(top.spleen[(spleen.blood.valid == 0 & spleen.non.valid == 0),][1001:1251,])
sum(spleen.valid)/length(spleen.valid)
sum(spleen.blood.valid)/length(spleen.blood.valid)
table(top.spleen[spleen.valid == 1,"Species"])
table(top.spleen[spleen.valid == 1,"GPL"])
gpl.name[match(names(table(top.spleen[spleen.valid == 1,"GPL"])), gpl.name$gpl),]
@
<<muscleQuantification>>=
muscle.meta<-cbind(tissueDistance$"skeletal muscle",
              GEO.metadata.matrix[
                match(rownames(tissueDistance$"skeletal muscle"),
                GEO.metadata.matrix$GSM)
              ,])
sum(muscle.meta$pvalue < 0.02)
sum(muscle.meta$distance < 0.2)

top.muscle<-muscle.meta[muscle.meta$pvalue < 0.02,]
top.muscle<-muscle.meta[muscle.meta$distance < 0.15,]
top.muscle<-muscle.meta[1:7500,]
non.muscle.terms<-c("kidney", "liver", "adipose")
# N.B. could use GSE12189 as a positive/negative control for muscle
muscle.valid<-(rowSums(
  sapply(muscle.terms, function(x){
    grepl(x, apply(
            top.muscle[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
muscle.non.valid<-(rowSums(
  sapply(non.muscle.terms, function(x){
    grepl(x, apply(
            top.muscle[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)



muscle.valid[top.muscle$GSE %in% muscle.GSE]<-1
non.muscle<-top.muscle[muscle.valid == 0,]
View(top.muscle[(muscle.valid == 0 & muscle.non.valid == 0),])
View(non.muscle)

sum(muscle.valid)/length(muscle.valid)
sum(muscle.valid[1:3000])/length(muscle.valid[1:3000])

table(top.muscle[muscle.valid == 1,"Species"])
table(top.muscle[muscle.valid == 1,"GPL"])
gpl.name[match(names(table(top.muscle[muscle.valid == 1,"GPL"])), gpl.name$gpl),]
@
<<kidneyQuantification>>=
kidney.meta<-cbind(tissueDistance$"kidney",
              GEO.metadata.matrix[
                match(rownames(tissueDistance$"kidney"),
                GEO.metadata.matrix$GSM)
              ,])
sum(kidney.meta$pvalue < 0.02)
sum(kidney.meta$distance < 0.2)

top.kidney<-kidney.meta[kidney.meta$pvalue < 0.02,]
top.kidney<-kidney.meta[kidney.meta$distance < 0.15,]
top.kidney<-kidney.meta[1:7500,]
non.kidney.terms<-c("liver", "adipose", "spleen", "jejunal", "gallbladder", "lung",
                    "hepatocellular", "duodenum", "intestine", "colorectal", "colon",
                    "intestinal", "hepatocytes", "jejunum", "ileum", "heart", "aorta",
                    "Gastric")
# N.B. could use GSE12189 as a positive/negative control for kidney
kidney.valid<-(rowSums(
  sapply(kidney.terms, function(x){
    grepl(x, apply(
            top.kidney[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
kidney.non.valid<-(rowSums(
  sapply(non.kidney.terms, function(x){
    grepl(x, apply(
            top.kidney[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)



kidney.valid[top.kidney$GSE %in% kidney.GSE]<-1
non.kidney<-top.kidney[kidney.valid == 0,]
View(top.kidney[(kidney.valid == 0 & kidney.non.valid == 0),])
View(non.kidney)

sum(kidney.valid)/length(kidney.valid)
sum(kidney.valid[1:3000])/length(kidney.valid[1:3000])


table(top.kidney[kidney.valid == 1,"Species"])
table(top.kidney[kidney.valid == 1,"GPL"])
gpl.name[match(names(table(top.kidney[kidney.valid == 1,"GPL"])), gpl.name$gpl),]
@
<<combine>>=
tissues.valid<-list(brain.valid, kidney.valid, liver.valid, lung.valid, muscle.valid, spleen.valid)
tissues.meta<-list(brain.meta, kidney.meta, liver.meta, lung.meta, muscle.meta, spleen.meta)
names(tissues.valid)<-tissues
names(tissues.meta)<-tissues
@