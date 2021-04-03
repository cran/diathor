#' Calculates the ILM Index (ILM)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Leclercq, L., & Maquet, B. (1987). Deux nouveaux indices diatomique et de qualité chimique des eaux courantes. Comparaison avec différents indices existants. Cahier de Biology Marine, 28, 303-310.
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
#' ilmResults <- diat_ilm(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_ilm


###### ---------- FUNCTION FOR ILM INDEX (Leclercq & Maq. 1988)---------- ########
### INPUT: resultLoad Data
### OUTPUTS: dataframe with ILM index per sample
diat_ilm <- function(resultLoad){

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
  #ilmDB <- read.csv("../Indices/ilm.csv") #uses the external csv file
  ilmDB <- diathor::ilm

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)


  #exact matches species in input data to acronym from index
  # taxaIn$ilm_v <- ilmDB$ilm_v[match(taxaIn$acronym, trimws(ilmDB$acronym))]
  # taxaIn$ilm_s <- ilmDB$ilm_s[match(taxaIn$acronym, trimws(ilmDB$acronym))]

  # #the ones still not found (NA), try against fullspecies
  taxaIn$ilm_v <- NA
  taxaIn$ilm_s <- NA
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$ilm_s[i]) | is.na(taxaIn$ilm_v[i])){
      taxaIn$ilm_v[i] <- ilmDB$ilm_v[match(trimws(rownames(taxaIn[i,])), trimws(ilmDB$fullspecies))]
      taxaIn$ilm_s[i] <- ilmDB$ilm_s[match(trimws(rownames(taxaIn[i,])), trimws(ilmDB$fullspecies))]
    }
  }


  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------ILM INDEX START --------#############
  print("Calculating ILM index")
  #creates results dataframe
  ilm.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(ilm.results) <- c("ILM", "ILM20", "Precision")
  #finds the column
  ilm_s <- (taxaIn[,"ilm_s"])
  ilm_v <- (taxaIn[,"ilm_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    ILMtaxaused <- (length(which(ilm_s * taxaIn[,sampleNumber] > 0))*100 / length(ilm_s))
    #print(paste("Total taxa input for ILM:", length(ilm_s)))
    #print(paste("Total taxa used for ILM:", sum(!is.na(taxaIn$ilm_s ))))
    #remove the NA
    ilm_s[is.na(ilm_s)] = 0
    ilm_v[is.na(ilm_v)] = 0
    ILM <- sum((taxaIn[,sampleNumber]*as.double(ilm_s)*as.double(ilm_v)))/sum(taxaIn[,sampleNumber]*as.double(ilm_v)) #raw value
    ILM20 <- (4.75*ILM)-3.75
    ilm.results[sampleNumber, ] <- c(ILM, ILM20,ILMtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------ILM INDEX: END--------############
  #PRECISION
  resultsPath <- resultLoad[[4]]
  precisionmatrix <- read.csv(paste(resultsPath,"\\Precision.csv", sep=""))
  precisionmatrix <- cbind(precisionmatrix, ilm.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="ilm.results$Precision"] <- "ILM"
  write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$ilm_s > 0)]
  inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv", sep=""))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "ILM")
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

  rownames(ilm.results) <- resultLoad[[3]]
  return(ilm.results)
}

