#' Calculates the SPEAR(herbicides) Index (SPEAR)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Wood, R. J., Mitrovic, S. M., Lim, R. P., Warne, M. S. J., Dunlop, J., & Kefford, B. J. (2019). Benthic diatoms as indicators of herbicide toxicity in rivers–A new SPEcies At Risk (SPEARherbicides) index. Ecological Indicators, 99, 203-213.
#' }
#'
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
#' }
#' @examples
#' \donttest{
#' # Example using sample data included in the package (sampleData):
#' data("diat_sampleData")
#' # First, the diat_loadData() function has to be called to read the data
#' # The data will be stored into a list (loadedData)
#' # And an output folder will be selected through a dialog box if resultsPath is empty
#' # In the example, a temporary directory will be used in resultsPath
#' df <- diat_loadData(diat_sampleData, resultsPath = tempdir())
#' spearResults <- diat_spear(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_spear

###### ---------- FUNCTION FOR SPEAR INDEX (Wood et al. 2019) ---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with SPEAR index per sample
diat_spear <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaInRA <- resultLoad[[1]]

  #Loads the species list specific for this index
  #spearDB <- read.csv("../Indices/spear.csv") #uses the external csv file
  spearDB <- diathor::spear

  #creates a species column with the rownames to fit in the script
  taxaInRA$species <- row.names(taxaInRA)

  #remove the "-"
  for (i in 1:nrow(spearDB)){
    if (spearDB$spear_v[i]=="-"){
      spearDB$spear_v[i] <- ""
    }
  }

  #exact matches species in input data to acronym from index
  taxaInRA$spear_v <- as.integer(spearDB$spear_v[match(taxaInRA$acronym, trimws(spearDB$acronym))])



  #the ones still not found (NA), try against fullspecies
  for (i in 1:nrow(taxaInRA)) {
    if (is.na(taxaInRA$spear_v[i])){
      taxaInRA$spear_v[i] <- spearDB$spear_v[match(trimws(rownames(taxaInRA[i,])), trimws(spearDB$fullspecies))]
    }
  }

  #removes NA from taxaInRA
  # taxaInRA[is.na(taxaInRA)] <- 0

  #gets the column named "acronym", everything before that is a sample
  lastcol <- which(colnames(taxaInRA)=="acronym")

  #######--------SPEAR INDEX START --------#############
  print("Calculating SPEAR index")
  #creates results dataframe
  spear.results <- data.frame(matrix(ncol = 2, nrow = (lastcol-1)))
  colnames(spear.results) <- c("SPEAR", "Precision")
  #finds the column
  spear_v <- (taxaInRA[,"spear_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    SPEARtaxaused <- (length(which(!is.na(spear_v))) * 100 / length(spear_v))
    SPEAR <- sum((log10(taxaInRA[,sampleNumber]+1)*as.double(spear_v)), na.rm = TRUE)/sum(log10(taxaInRA[,sampleNumber]+1), na.rm = TRUE) #raw value
    SPEAR <- SPEAR * 100
    spear.results[sampleNumber, ] <- c(SPEAR, SPEARtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------SPEAR INDEX: END--------############
  #PRECISION
  resultsPath <- resultLoad[[4]]
  precisionmatrix <- read.csv(paste(resultsPath,"\\Precision.csv", sep=""))
  precisionmatrix <- cbind(precisionmatrix, spear.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="spear.results$Precision"] <- "SPEAR"
  write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))
  #END PRECISION



  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaInRA$species[which(!is.na(taxaInRA$spear_v))]
  inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv", sep=""))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "SPEAR")
  #creates a new data matrix to append the existing Taxa Included file
  newinclusionmatrix <- as.data.frame(matrix(nrow=max(length(taxaIncluded), nrow(inclusionmatrix)), ncol=ncol(inclusionmatrix)+1))
  for (i in 1:ncol(inclusionmatrix)){
    newinclusionmatrix[1:nrow(inclusionmatrix),i] <- as.character(inclusionmatrix[1:nrow(inclusionmatrix),i])
  }
  if (nrow(newinclusionmatrix) > length(taxaIncluded)){
    newinclusionmatrix[1:length(taxaIncluded), ncol(newinclusionmatrix)] <- taxaIncluded
  } else {
    newinclusionmatrix[1:nrow(newinclusionmatrix), ncol(newinclusionmatrix)] <- taxaIncluded
  }
  inclusionmatrix <- newinclusionmatrix
  colnames(inclusionmatrix) <- colnamesInclusionMatrix
  inclusionmatrix <- inclusionmatrix[-(1:which(colnames(inclusionmatrix)=="Eco.Morpho")-1)]
  write.csv(inclusionmatrix, paste(resultsPath,"\\Taxa included.csv", sep=""))
  #END TAXA INCLUSION
  #EXCLUDED TAXA
  taxaExcluded <- taxaInRA[!('%in%'(taxaInRA$species,taxaIncluded)),"species"]
  exclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  #creates a new data matrix to append the existing Taxa Included file
  newexclusionmatrix <- as.data.frame(matrix(nrow=max(length(taxaExcluded), nrow(exclusionmatrix)), ncol=ncol(exclusionmatrix)+1))
  for (i in 1:ncol(exclusionmatrix)){
    newexclusionmatrix[1:nrow(exclusionmatrix),i] <- as.character(exclusionmatrix[1:nrow(exclusionmatrix),i])
  }
  if (nrow(newexclusionmatrix) > length(taxaExcluded)){
    newexclusionmatrix[1:length(taxaExcluded), ncol(newexclusionmatrix)] <- taxaExcluded
  } else {
    newexclusionmatrix[1:nrow(newexclusionmatrix), ncol(newexclusionmatrix)] <- taxaExcluded
  }
  exclusionmatrix <- newexclusionmatrix
  colnames(exclusionmatrix) <- colnamesInclusionMatrix
  exclusionmatrix <- exclusionmatrix[-(1:which(colnames(exclusionmatrix)=="Eco.Morpho")-1)]
  write.csv(exclusionmatrix, paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  #END EXCLUDED TAXA

  rownames(spear.results) <- resultLoad[[3]]
  return(spear.results)
}



