#' Calculates the Pampean Diatom Index (IDP)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Gómez, N., & Licursi, M. (2001). The Pampean Diatom Index (IDP) for assessment of rivers and streams in Argentina. Aquatic Ecology, 35(2), 173-181.
#' }
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
#' idpResults <- diat_idp(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_idp


###### ---------- FUNCTION FOR IDP INDEX (Pampean Index - Gomez & Licursi)  ---------- ########
#### IN THIS SECTION WE CALCULATE IDP INDEX (Pampean Index - Gomez & Licursi)
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with IDP index per sample
diat_idp <- function(resultLoad){

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
  #idpDB <- read.csv("../Indices/idp.csv") #uses the external csv file
  idpDB <- diathor::idp

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  #if acronyms exist, use them, its more precise
  #if there is an acronym column, it removes it and stores it for later
  #exact matches species in input data to acronym from index
  # taxaIn$idp_v <- idpDB$idp_v[match(trimws(taxaIn$acronym), trimws(idpDB$acronym))]

  taxaIn$idp_v <- NA
  #the ones still not found (NA), try against fullspecies
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$idp_v[i])){
      taxaIn$idp_v[i] <- idpDB$idp_v[match(trimws(rownames(taxaIn[i,])), trimws(idpDB$fullspecies))]
    }
  }


  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------IDP INDEX START --------#############
  print("Calculating IDP index")
  idp.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(idp.results) <- c("IDP", "IDP20", "Precision")
  #finds the column
  idp_v <- (taxaIn[,"idp_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    IDPtaxaused <- (length(which(idp_v * taxaIn[,sampleNumber] > 0))*100 / length(idp_v))

    #print(paste("Total taxa used for IDP:", sum(!is.na(taxaIn$idp_v ))))
    #remove the NA
    idp_v[is.na(idp_v)] = 0
    IDP <- sum((taxaIn[,sampleNumber]*as.double(idp_v)))/sum(taxaIn[which(idp_v > 0),sampleNumber]) #raw value
    IDP20 <- 20-(4.75*IDP)
    idp.results[sampleNumber, ] <- c(IDP, IDP20,IDPtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #outputs: IDP, IDP20
  #close progressbar
  close(pb)
  #######--------IDP INDEX: END--------############
  #PRECISION
  resultsPath <- resultLoad[[4]]

  #precisionmatrix <- read.csv(paste(resultsPath,"\\Precision.csv", sep=""))
  precisionmatrix <- read.csv(file.path(resultsPath, "Precision.csv"))
  precisionmatrix <- cbind(precisionmatrix, idp.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="idp.results$Precision"] <- "IDP"
  #write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))
  write.csv(precisionmatrix, file.path(resultsPath, "Precision.csv"))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$idp_v > 0)]
  #inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv", sep=""))
  inclusionmatrix <- read.csv(file.path(resultsPath, "Taxa included.csv"))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "IDP")
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
  taxaExcluded <- taxaIn[!('%in%'(taxaIn$species,taxaIncluded)),"species"]
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

  rownames(idp.results) <- resultLoad[[3]]
  return(idp.results)
}

