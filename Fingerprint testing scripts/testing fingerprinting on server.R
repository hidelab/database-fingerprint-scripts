# Script for speed testing the fingerprint
source("/home/galtschu2/Documents/Databases/gabriel functions.R")
# Profile second time around as objects are loaded into the global environment the first time

Rprof(tmp <- tempfile())
test<-geo2fingerprint("GSM27832", "KEGG and Wikipathways and static", enrichmentMethod = "SCG")
Rprof()
summaryRprof(tmp)
unlink(tmp)

Rprof(tmp <- tempfile())
test1<-geo2fingerprint("GSM27832", "KEGG and Wikipathways and static", enrichmentMethod = "SCE")
Rprof()
summaryRprof(tmp)
unlink(tmp)

# 19.38s vs 5.38s - significant improvement
# compare 16.76 for single.chip.GSEA to 2.68 for single.chip.enrichment within the scripts
# most of the time lag function in single chip enrichemnt (2.68) and customCDFAnn (0.84)
# individual commands apply, mean, FUN and dir take up the most time
# dir is a command that retrieves the full list of GEO files in a directory
# if the GEO file is not found then it is downloaded
# Not necessary repeatedly do this, better to load the list into the global environment
# However, this means that this will have to be updated as well be the GEO loading script

# updated so that the GEO directory data is only loaded once into the Global Environment

Rprof(tmp <- tempfile())
test1<-geo2fingerprint("GSM27832", "KEGG and Wikipathways and static", enrichmentMethod = "SCE")
Rprof()
summaryRprof(tmp)
unlink(tmp)

# great, now down to 3.78 seconds
# also of note is that SCE returns a matrix while SCG returns a dataframe
# matrices may be more efficient

# updated script with alternative path to new chipframe, which includes Illumina arrays
# test Illumina array fingerprint
# Fingerprinting Illumina worked!

Rprof(tmp <- tempfile())
test2<-geo2fingerprint("GSM365199", "KEGG and Wikipathways and static", enrichmentMethod = "SCE")
Rprof()
summaryRprof(tmp)
unlink(tmp)

Rprof(tmp <- tempfile())
test3<-geo2fingerprint("GSM365199", "KEGG and Wikipathways and static", enrichmentMethod = "SCG")
Rprof()
summaryRprof(tmp)
unlink(tmp)

