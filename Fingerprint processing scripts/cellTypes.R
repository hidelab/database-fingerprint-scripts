# Compiling cell types
# Author: Gabriel Altschuler
# Timestamp: 20110629
# Status: complete
# Script to compile a dataframe of the curated cell types

load(
 "/data/shared/Fingerprint/curatedCellTypes/curatedPluripotentArrays_20110510.RData"
    )
load(
 "/data/shared/Fingerprint/curatedCellTypes/curatedESArrays.RData"
    )
load(
 "/data/shared/Fingerprint/curatedCellTypes/curatediPSArrays.RData"
    )
load(
 "/data/shared/Fingerprint/curatedCellTypes/curatedFibroblastsArrays.RData"
    )
load(
 "/data/shared/Fingerprint/curatedCellTypes/curatedMelanomaSamples.RData"
    )
load(
 "/data/shared/Fingerprint/curatedCellTypes/putative.epithelial.RData"
    )
load(
 "/data/shared/Fingerprint/curatedCellTypes/putative.mesenchymal.stem.RData"
    )
load(
 "/data/shared/Fingerprint/curatedCellTypes/putative.neural.crest.RData"
    )

# don't include epithelial or fibroblasts for now

cellTypes<-data.frame(GSM = c(
  pluripotents.frame$GSM,
  ES.frame$GSM,
  iPS.frame$GSM,
  melanoma.frame$gsm,
  putative.mesenchymal.stem$GSM,
  putative.neural.crest$GSM),
                      Type = c(
  rep("pluripotent", length(pluripotents.frame$GSM)),
  rep("ES", length(ES.frame$GSM)),
  rep("iPS", length(iPS.frame$GSM)),
  rep("melanoma", length(melanoma.frame$gsm)),
  rep("mesenchymal", length(putative.mesenchymal.stem$GSM)),
  rep("neuralCrest", length(putative.neural.crest$GSM))
  ), stringsAsFactors = FALSE)

save(cellTypes, file = "/data/shared/Fingerprint/curatedCellTypes/cellTypes_collection.RData")