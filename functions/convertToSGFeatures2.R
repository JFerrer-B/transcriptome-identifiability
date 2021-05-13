convertToSGFeatures2 <- function (x, coerce = FALSE, merge = FALSE) 
{
  if (!is(x, "TxFeatures")) {
    stop("x must be a TxFeatures object")
  }
  if (length(x) == 0) {
    return(SGFeatures())
  }
  if (coerce) {
    features <- granges(x)
    mcols(features)$type <- as.character(type(x))
    splice5p <- mcols(features)$type %in% c("I", "L")
    splice3p <- mcols(features)$type %in% c("I", "F")
    splice5p[mcols(features)$type == "J"] <- NA
    splice3p[mcols(features)$type == "J"] <- NA
    mcols(features)$type[mcols(features)$type != "J"] <- "E"
    mcols(features)$splice5p <- splice5p
    mcols(features)$splice3p <- splice3p
    mcols(features)$txName <- txName(x)
    mcols(features)$geneName <- geneName(x)
  }
  else {
    features <- processFeatures2(x, merge = FALSE)
  }
  features <- SGSeq:::addFeatureID(features)
  features <- SGSeq:::addGeneID(features)
  features <- SGFeatures(features)
  if (!coerce) {
    features <- annotate(features, x)
  }
  return(features)
}