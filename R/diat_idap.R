#' Calculates the Indice Diatomique Artois-Picardie (IDAP)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Prygiel, J., & Coste, M. (1993). The assessment of water quality in the Artois-Picardie water basin (France) by the use of diatom indices. Hydrobiologia, 269(1), 343-349.
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
#' idapResults <- diat_idap(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_idap



###### ---------- FUNCTION FOR IDAP INDEX (Indice Diatomique Artois-Picardie) ---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with IDAP index per sample
diat_idap <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]] #indices use raw matrix

  ### START NEW CORRECTIONS
  #Loads the species list specific for this index
  #idapDB <- read.csv("../Indices/idap.csv") #uses the external csv file
  idapDB <- diathor::idap

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  # #exact matches species in input data to acronym from index
  # taxaIn$idap_v <- idapDB$idap_v[match(taxaIn$acronym, trimws(idapDB$acronym))]
  # taxaIn$idap_s <- idapDB$idap_s[match(taxaIn$acronym, trimws(idapDB$acronym))]

  # #the ones still not found (NA), try against fullspecies
  taxaIn$idap_v <- NA
  taxaIn$idap_s <- NA
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$idap_s[i]) | is.na(taxaIn$idap_v[i])){
      taxaIn$idap_v[i] <- idapDB$idap_v[match(trimws(rownames(taxaIn[i,])), trimws(idapDB$fullspecies))]
      taxaIn$idap_s[i] <- idapDB$idap_s[match(trimws(rownames(taxaIn[i,])), trimws(idapDB$fullspecies))]
    }
  }


  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------IDAP INDEX START --------#############
  print("Calculating IDAP index")
  #creates results dataframe
  idap.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(idap.results) <- c("IDAP", "IDAP20", "Precision")
  #finds the column
  idap_s <- (taxaIn[,"idap_s"])
  idap_v <- (taxaIn[,"idap_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    IDAPtaxaused <- (length(which(idap_s * taxaIn[,sampleNumber] > 0))*100 / length(idap_s))
    #remove the NA
    idap_s[is.na(idap_s)] = 0
    idap_v[is.na(idap_v)] = 0
    IDAP <- sum((taxaIn[,sampleNumber]*as.double(idap_s)*as.double(idap_v)))/sum(taxaIn[,sampleNumber]*as.double(idap_v)) #raw value
    IDAP20 <- (4.75*IDAP)-3.75
    idap.results[sampleNumber, ] <- c(IDAP, IDAP20,IDAPtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------IDAP INDEX: END--------############
  #PRECISION
  resultsPath <- resultLoad[[4]]
  #precisionmatrix <- read.csv(paste(resultsPath,"\\Precision.csv", sep=""))
  precisionmatrix <- read.csv(file.path(resultsPath, "Precision.csv"))
  precisionmatrix <- cbind(precisionmatrix, idap.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="idap.results$Precision"] <- "IDAP"
  #write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))
  write.csv(precisionmatrix, file.path(resultsPath, "Precision.csv"))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$new_species[which(taxaIn$idap_s > 0)]
  #inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv"))
  inclusionmatrix <- read.csv(file.path(resultsPath, "Taxa included.csv"))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "IDAP")
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
  #write.csv(inclusionmatrix, paste(resultsPath,"\\Taxa included.csv", sep=""))
  write.csv(inclusionmatrix, file.path(resultsPath,"Taxa included.csv"))
  #END TAXA INCLUSION
  #EXCLUDED TAXA
  taxaExcluded <- taxaIn[!('%in%'(taxaIn$new_species,taxaIncluded)),"new_species"]
  #exclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  exclusionmatrix <- read.csv(file.path(resultsPath, "Taxa excluded.csv"))
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
  #write.csv(exclusionmatrix, paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  write.csv(exclusionmatrix, file.path(resultsPath,"Taxa excluded.csv"))
  #END EXCLUDED TAXA

  rownames(idap.results) <- resultLoad[[3]]
  return(idap.results)
}



