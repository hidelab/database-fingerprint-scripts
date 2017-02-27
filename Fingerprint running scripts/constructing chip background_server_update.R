setwd("/home/galtschu2/Documents/Projects/Fingerprinting")
source("/home/galtschu2/Documents/Databases/gabriel functions.R")
# stopped saving workspace as files very large

GPL570.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU133Plus2_Hs_ENTREZG_mapping.txt")
GPL1261.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/Mouse4302_Mm_ENTREZG_mapping.txt")
GPL339.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/MOE430A_Mm_ENTREZG_mapping.txt")
GPL96.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU133A_Hs_ENTREZG_mapping.txt")
GPL97.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU133B_Hs_ENTREZG_mapping.txt")
GPL81.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/MGU74Av2_Mm_ENTREZG_mapping.txt")
GPL8321.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/Mouse430A2_Mm_ENTREZG_mapping.txt")
GPL8300.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU95Av2_Hs_ENTREZG_mapping.txt")
GPL571.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/HGU133A2_Hs_ENTREZG_mapping.txt")
GPL340.data<-read.delim(file = "/home/galtschu2/Documents/Databases/customCDFmappings/MOE430B_Mm_ENTREZG_mapping.txt")
load("/home/galtschu2/Documents/Databases/customCDFmappings/GPL2986.ann.RData")
chiplist<-c("GPL570", "GPL1261", "GPL339", "GPL96", "GPL97", "GPL81", "GPL8321", "GPL8300", "GPL571", "GPL340", "GPL2986")

annotation<-list(
GPL570.ann<-annotation(GPL570.data),
GPL1261.ann<-annotation(GPL1261.data),
GPL339.ann<-annotation(GPL339.data),
GPL96.ann<-annotation(GPL96.data),
GPL97.ann<-annotation(GPL97.data),
GPL81.ann<-annotation(GPL81.data),
GPL8321.ann<-annotation(GPL8321.data),
GPL8300.ann<-annotation(GPL8300.data),
GPL571.ann<-annotation(GPL571.data),
GPL340.ann<-annotation(GPL340.data),
GPL2986.ann<-GPL2986.ann
)
names(annotation) <- chiplist


chipframe = vector("list", length = 0)
for (i in 1:length(chiplist)){
	chipframe <- append(chipframe, list(list(ann = annotation[[chiplist[i]]])))
	}
names(chipframe)<-chiplist

save(chipframe, file = "data/chipframe_update.RData")

