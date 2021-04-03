#' Calculates the Swiss Diatom Index (IDCH)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Hürlimann J., Niederhauser P. 2007: Méthodes d’analyse et d’appréciation des cours d’eau. Diatomées Niveau R (région). État de l’environnement n° 0740. Office fédéral de l’environnement, Berne. 132 p
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
#' idchResults <- diat_idch(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8

#' @export diat_idch


###### ---------- FUNCTION FOR IDCH INDEX (Swiss Diatom Index, Lecointe et al., 2003)  ---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with IDCH index per sample
diat_idch <- function(resultLoad){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #idchDB <- read.csv("../Indices/idch.csv") #uses the external csv file
  idchDB <- diathor::idch
  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  # #exact matches species in input data to acronym from index
  # taxaIn$idch_v <- idchDB$idch_v[match(taxaIn$acronym, trimws(idchDB$acronym))]
  # taxaIn$idch_s <- idchDB$idch_s[match(taxaIn$acronym, trimws(idchDB$acronym))]

  # #the ones still not found (NA), try against fullspecies
  taxaIn$idch_v <- NA
  taxaIn$idch_s <- NA
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$idch_s[i]) | is.na(taxaIn$idch_v[i])){
      taxaIn$idch_v[i] <- idchDB$idch_v[match(trimws(rownames(taxaIn[i,])), trimws(idchDB$fullspecies))]
      taxaIn$idch_s[i] <- idchDB$idch_s[match(trimws(rownames(taxaIn[i,])), trimws(idchDB$fullspecies))]
    }
  }

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------IDCH INDEX START --------#############
  print("Calculating IDCH index")
  #creates results dataframe
  idch.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(idch.results) <- c("IDCH", "IDCH20", "Precision")
  #finds the column
  idch_s <- (taxaIn[,"idch_s"])
  idch_v <- (taxaIn[,"idch_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    IDCHtaxaused <- (length(which(idch_s * taxaIn[,sampleNumber] > 0))*100 / length(idch_s))
    #remove the NA
    idch_s[is.na(idch_s)] = 0
    idch_v[is.na(idch_v)] = 0
    IDCH <- sum((taxaIn[,sampleNumber]*as.double(idch_s)*as.double(idch_v)))/sum(taxaIn[,sampleNumber]*as.double(idch_v)) #raw value
    IDCH20 <- 22.714-(2.714*IDCH)
    idch.results[sampleNumber, ] <- c(IDCH, IDCH20,IDCHtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------IDCH INDEX: END--------############
  #PRECISION
  resultsPath <- resultLoad[[4]]
  precisionmatrix <- read.csv(paste(resultsPath,"\\Precision.csv", sep=""))
  precisionmatrix <- cbind(precisionmatrix, idch.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="idch.results$Precision"] <- "IDCH"
  write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$idch_s > 0)]
  inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv", sep=""))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "IDCH")
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
  taxaExcluded <- taxaIn[!('%in%'(taxaIn$species,taxaIncluded)),"species"]
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

  rownames(idch.results) <- resultLoad[[3]]
  return(idch.results)

}

