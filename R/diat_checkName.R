#' Searches all the taxa database for the input name
#' @param taxaname the name of the taxa (genus, species, variety) to be checked against the internal DB
#' @param byword if byword = F (default), the input string will be searched without splitting words. If True, each word will be searched separately
#' @description
#' Searches all the taxa database for the input name, returns a list with the results
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @importFrom stringdist stringdist ain
#' @export diat_checkName
#'

diat_checkName <- function(taxaname, byword = F) {
  #load databases
  getDiatBarcode <- diat_getDiatBarcode()
  ecodata <- as.data.frame(getDiatBarcode[1]) #ecodata
  taxaList <- as.data.frame(getDiatBarcode[2]) #ecodata

  #result List
  resultList <- NA

  #if byword = T, separates the input taxaname by word and searches them independently
  if (byword == F){
    searchvectr <- sort(taxaList$species[str_detect(taxaList$species, taxaname)])
    resultList <- searchvectr
  } else {
    for (i in 1:lengths(strsplit(taxaname, " "))){
      searchvectr <- taxaList$species[str_detect(taxaList$species, word(taxaname,i))]
      resultList <- c(resultList, searchvectr)
    }
  }
  return(sort(resultList))
}
