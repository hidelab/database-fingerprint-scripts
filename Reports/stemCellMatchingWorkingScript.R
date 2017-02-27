# now try to find 'interesting datasets
sum(head(typeDistance.meta[["ES"]]$pluripotent,1000))
# 730 of 1000 are pluripotent
temp1<-head(typeDistance.meta[["ES"]],1000)
temp2<-temp1[temp1$pluripotent == 0,]
temp3<-temp2[!(grepl("embryo", apply(
            temp2[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)),]

View(temp3[!(grepl("fetal", apply(
            temp3[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)),])


dim(temp2[(grepl("embryo", apply(
            temp2[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)),])

temp<-head(typeDistance.meta[["ES"]][typeDistance.meta[["ES"]]$pluripotent == 0,],1000)
View(temp)
View(head(temp[order(temp$GSEcount, decreasing = TRUE),],20))

<<matchingArrays>>=
top.pluripotent<-pluripotent.meta[pluripotent.meta$distance < 0.1,]
top.ES<-ES.meta[ES.meta$distance < 0.1,]
top.iPS<-iPS.meta[iPS.meta$distance < 0.1,]

                

pluripotent.valid<-(rowSums(
  sapply(pluripotent.terms, function(x){
    grepl(x, apply(
            top.pluripotent[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)


ES.valid<-(rowSums(
  sapply(ES.terms, function(x){
    grepl(x, apply(
            top.ES[, c("Title", "Source", "Characteristics")], 1, function(x){
                      paste(x, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
                
                

                
pluripotent.valid[top.pluripotent$GSE %in% pluripotent.GSE]<-1
non.pluripotent<-top.pluripotent[pluripotent.valid == 0,]
non.ES<-top.ES[ES.valid == 0,]

View(non.pluripotent)
View(non.ES)

EXPO dataset
GSE2109

hist(typeDistance.meta[["ES"]][GEO.metadata.matrix$GSM[GEO.metadata.matrix$GSE %in% "GSE2109"], "distance"], breaks = seq(0,1,1/sum(abs(typeConsensus[,"ES"]))))

<<test>>=
hist(pluripotent.distance[GEO.metadata.matrix$GSM[GEO.metadata.matrix$GSE %in% "GSE19274"], "distance"], breaks = seq(0,1,1/sum(abs(pluripotent.consensus))))
hist(ES.distance[GEO.metadata.matrix$GSM[GEO.metadata.matrix$GSE %in% "GSE19274"], "distance"], breaks = seq(0,1,1/sum(abs(ES.consensus))))

GSE19274.ES<-ES.meta[ES.meta$GSE %in% "GSE19274",]

GSE19274.ES$tissue<-sapply(GSE19274.ES$Characteristics, function(x){unlist(strsplit(x, split = "\t"))[2]})
GSE19274.ES$risk<-sapply(GSE19274.ES$Characteristics, function(x){unlist(strsplit(x, split = "\t"))[6]})
GSE19274.ES$mycn<-sapply(GSE19274.ES$Characteristics, function(x){unlist(strsplit(x, split = "\t"))[3]})


library(beanplot)
beanplot(distance ~ risk, data = GSE19274.ES[grep("primary",GSE19274.ES$tissue),])
beanplot(distance ~ mycn, data = GSE19274.ES[grep("primary",GSE19274.ES$tissue),])
beanplot(distance ~ tissue, data = GSE19274.ES[grep("primary",GSE19274.ES$tissue),])




hist(pluripotent.distance[GEO.metadata.matrix$GSM[GEO.metadata.matrix$GSE %in% "GSE10843"], "distance"], breaks = seq(0,1,1/sum(abs(pluripotent.consensus))))
hist(pluripotent.distance[GEO.metadata.matrix$GSM[GEO.metadata.matrix$GSE %in% "GSE23806"], "distance"], breaks = seq(0,1,1/sum(abs(pluripotent.consensus))))


hist(pluripotent.distance[GEO.metadata.matrix$GSM[GEO.metadata.matrix$GSE %in% "GSE23806"], "distance"], breaks = seq(0,1,1/sum(abs(pluripotent.consensus))))


tail(GEO.metadata.matrix[(GEO.metadata.matrix$GSE == "GSE23806"),])

GSE12276.ES<-typeDistance.meta[["ES"]][typeDistance.meta[["ES"]]$GSE %in% "GSE12276",]
GSE12276.ES$survival<-as.numeric(gsub("Survival time \\(months\\): ", "", sapply(GSE12276.ES$Characteristics, function(x){unlist(strsplit(x, split = "\t"))[[2]]})))
GSE12276.ES$distGroup<-kmeans(GSE12276.ES$distance, 2)$cluster
fit<-survfit(Surv(GSE12276.ES$survival) ~ GSE12276.ES$distGroup)
hist(GSE12276.ES$distance)

GSE4271.ES<-typeDistance.meta[["ES"]][typeDistance.meta[["ES"]]$GSE %in% "GSE4271",]
GSE4271.ES$survival<-as.numeric(gsub(" Survival \\(weeks\\): ", "", sapply(GSE4271.ES$Characteristics, function(x){unlist(strsplit(x, split = ";"))[[4]]})))
GSE4271.ES$censor<-gsub(" Censor: ", "", sapply(GSE4271.ES$Characteristics, function(x){unlist(strsplit(x, split = ";"))[[5]]}))
GSE4271.ES$censor<-gsub("no", 1, GSE4271.ES$censor)
GSE4271.ES$censor<-gsub("yes", 0, GSE4271.ES$censor)
GSE4271.ES$censor<-as.numeric(GSE4271.ES$censor)
GSE4271.ES$distGroup<-kmeans(GSE4271.ES$distance, 2)$cluster
fit<-survfit(Surv(GSE4271.ES$survival, GSE4271.ES$censor) ~ GSE4271.ES$distGroup)
plot(fit)
hist(GSE4271.ES$distance)
plot(GSE4271.ES$distance, GSE4271.ES$survival)

# this one works okay
GSE4412.ES<-typeDistance.meta[["ES"]][typeDistance.meta[["ES"]]$GSE %in% "GSE4412",]
GSE4412.ES$survival<-as.numeric(gsub("survival time: ", "", sapply(GSE4412.ES$Characteristics, function(x){unlist(strsplit(x, split = ";\t"))[[8]]})))
GSE4412.ES$censor<-gsub("living: ", "", sapply(GSE4412.ES$Characteristics, function(x){unlist(strsplit(x, split = ";\t"))[[7]]}))
GSE4412.ES$censor<-gsub("DECEASED", 1, GSE4412.ES$censor)
GSE4412.ES$censor<-gsub("ALIVE", 0, GSE4412.ES$censor)
GSE4412.ES$censor<-as.numeric(GSE4412.ES$censor)
GSE4412.ES$distGroup<-kmeans(GSE4412.ES$distance, 2)$cluster
fit<-survfit(Surv(GSE4412.ES$survival, GSE4412.ES$censor) ~ GSE4412.ES$distGroup)
plot(fit)
hist(GSE4412.ES$distance)
plot(GSE4412.ES$distance, GSE4412.ES$survival)

gpl.description<-dbGetQuery(con, "select gsm, description from gsm")
GSE2109.ES<-typeDistance.meta[["ES"]][typeDistance.meta[["ES"]]$GSE %in% "GSE2109",]
GSE2109.ES$description<-gpl.description[match(rownames(GSE2109.ES), gpl.description$gsm),2]
GSE2109.ES<-GSE2109.ES[(grepl("Histology", GSE2109.ES$description)),]
GSE2109.ES$histology<-gsub(": ", "", sapply(GSE2109.ES$description, function(x){unlist(strsplit(x, split = "Histology"))[[2]]}))
GSE2109.ES$grade<-gsub(": ", "", sapply(GSE2109.ES$description, function(x){unlist(strsplit(x, split = "Clinical Grade: "))[[2]]}))
GSE2109.ES$tissue<-sapply(GSE2109.ES$Title, function(x){unlist(strsplit(x, split = " - "))[[1]]})

# try gene chip oncology database

large<-names(table(GSE2109.ES$histology)[table(GSE2109.ES$histology) > 20])
large.tissue<-names(table(GSE2109.ES$tissue)[table(GSE2109.ES$tissue) > 20])
GSE2109.ES$large<-GSE2109.ES$histology %in% large
GSE2109.ES$largeTissue<-GSE2109.ES$tissue %in% large.tissue
boxplot()
boxplot(distance ~ histology, data = GSE2109.ES, subset = large == TRUE, las = 2)
boxplot(distance ~ tissue, data = GSE2109.ES, subset = largeTissue == TRUE, las = 2, main = "expO cancer data - dist to ES consensus")

# ovarian cancer

GSE14407.ES<-typeDistance.meta[["ES"]][typeDistance.meta[["ES"]]$GSE %in% "GSE14407",]
GSE14407.ES$class<-sapply(GSE14407.ES$Characteristics, function(x){unlist(strsplit(x, split = "specimen: "))[[2]]})
boxplot(distance ~ class, data = GSE14407.ES, las = 2, main = "Ovarian cancer data - dist to ES consensus")

par(mar = c(15,5,5,5))


<<defineTerms2>>=
ovarianCancer.terms<-c("ovarian cancer")
breastCancer.terms<-c("breast cancer")
ewing.terms<-c("ewing")
lungCancer.terms<-c("lung cancer")
leukemia.terms<-c("leukemia")
endometrialCancer.terms<-c("Endometrial ancer")
cancer.terms<-list(
  ovarianCancer.terms,
  breastCancer.terms,
  ewing.terms,
  lungCancer.terms,
  leukemia.terms,
  endometrialCancer.terms
  )
names(cancer.terms)<-c("ovarianCancer", "breastCancer",
                       "ewing", "lungCancer", "leukemia",
                       "endometrialCancer")

ES.cancer.ann<-typeDistance.meta[["ES"]]

con <- dbConnect(SQLite(),
 "/Users/GabrielAltschuler/Documents/Databases/GEOmetadb/GEOmetadb.sqlite"
   )

gse.title<-dbGetQuery(con, "select gse, title from gse")
dbDisconnect(con)
ES.cancer.ann$GSEtitle<-gse.title[match(ES.cancer.ann$GSE, gse.title$gse),2]

for (i in names(cancer.terms)){
print(i)
  ES.cancer.ann[,i]<-(rowSums(
  sapply(cancer.terms[[i]], function(x){
    grepl(x, apply(
            ES.cancer.ann[, c("Title", "Source", "Characteristics", "GSEtitle")],
            1, function(y){
                      paste(y, collapse = " ")
                      }
                  ), ignore.case = TRUE)
      })
  ) > 0)
}

ES.cancer.ann$type<-NA
for (i in names(cancer.terms)){
  ES.cancer.ann$type[ES.cancer.ann[,i]]<-i
}



