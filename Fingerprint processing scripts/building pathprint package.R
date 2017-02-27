# building pathprint R package
scripts<-c("exprs2fingerprint",
          "consensusDistance",
          "consensusFingerprint",
          "customCDFAnn",
          "diffPathways",
          "single.chip.enrichment",
          "thresholdFingerprint"
          )
data<-c("chipframe",
        "GEO.fingerprint.matrix",
        "GEO.metadata.matrix",
        "platform.thresholds",
        "pluripotents.frame",
        "genesets"
        )
genesets<-c("pathprint.v0.3.Ce.gs",
            "pathprint.v0.3.Dm.gs",
            "pathprint.v0.3.Dr.gs",
            "pathprint.v0.3.Hs.gs",
            "pathprint.v0.3.Mm.gs",
            "pathprint.v0.3.Rn.gs"
            )


setwd("/Users/GabrielAltschuler/Dropbox/fingerprint")
sapply(scripts, function(x){source(paste("scripts/", x, ".R", sep = ""))})
sapply(data, function(x){load(paste("data/", x, ".RData", sep = ""), .GlobalEnv)})
sapply(genesets, function(x){load(paste("data/", x, sep = ""), .GlobalEnv)})

setwd("/Users/GabrielAltschuler/Documents/Projects/Fingerprinting/R_package/testArea/testPackage")
# Overwrite previous package directory if present
package.skeleton(list = c(scripts, data, genesets), name = "pathprint.v0.3.beta3", force = TRUE)
setwd("pathprint.v0.3.beta3/man")
man<-dir()
man<-man[-grep("package.Rd", man)]
for (i in 1:length(man)){
  code<-paste("sh ../../shelltest/mansimple.sh", man[i], sep = " ")
  system(code)
}
setwd("../..")
# system("R CMD check pathprint.v0.3.beta3")
system("R CMD build pathprint.v0.3.beta3")


# N.B. installed on server using
# install.packages("/home/galtschu2/Documents/Projects/Fingerprinting/R_package/pathprint.v0.3.beta3_1.0.tar.gz", lib = .libPaths()[3], repos = NULL, type = "source")
