#' Calculates the Lobo Index (LOBO)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Lobo, E. A., Callegaro, V. L. M., & Bender, E. P. (2002). Utilização de algas diatomáceas epilíticas como indicadoras da qualidade da água em rios e arroios da Região Hidrográfica do Guaíba, RS, Brasil. Edunisc.
#' }
#' \itemize{
#' \item Lobo, E. A., Bes, D., Tudesque, L., & Ector, L. (2004). Water quality assessment of the Pardinho River, RS, Brazil, using epilithic diatom assemblages and faecal coliforms as biological indicators. Vie et Milieu, 54(2-3), 115-126.
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
#' loboResults <- diat_lobo(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_lobo

###### ---------- FUNCTION FOR LOBO INDEX (Lobo et al. 2002)---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with LOBO index per sample
diat_lobo <- function(resultLoad){

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
  #loboDB <- read.csv("../Indices/lobo.csv") #uses the external csv file
  loboDB <- diathor::lobo


  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)


  # #exact matches species in input data to acronym from index
  # taxaIn$lobo_v <- loboDB$lobo_v[match(taxaIn$acronym, trimws(loboDB$acronym))]
  # taxaIn$lobo_s <- loboDB$lobo_s[match(taxaIn$acronym, trimws(loboDB$acronym))]

  # #the ones still not found (NA), try against fullspecies
  taxaIn$lobo_v <- NA
  taxaIn$lobo_s <- NA
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$lobo_s[i]) | is.na(taxaIn$lobo_v[i])){
      taxaIn$lobo_v[i] <- loboDB$lobo_v[match(trimws(rownames(taxaIn[i,])), trimws(loboDB$fullspecies))]
      taxaIn$lobo_s[i] <- loboDB$lobo_s[match(trimws(rownames(taxaIn[i,])), trimws(loboDB$fullspecies))]
    }
  }


  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------LOBO INDEX START --------#############
  print("Calculating LOBO index")
  #creates results dataframe
  lobo.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(lobo.results) <- c("LOBO", "LOBO20", "Precision")
  #finds the column
  lobo_s <- (taxaIn[,"lobo_s"])
  lobo_v <- (taxaIn[,"lobo_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    LOBOtaxaused <- (length(lobo_s) - sum(is.na(taxaIn$lobo_s )))*100 / length(lobo_s)
    LOBOtaxaused <- (length(which(lobo_s * taxaIn[,sampleNumber] > 0))*100 / length(lobo_s))
    #remove the NA
    lobo_s[is.na(lobo_s)] = 0
    lobo_v[is.na(lobo_v)] = 0
    LOBO <- sum((taxaIn[,sampleNumber]*as.double(lobo_s)*as.double(lobo_v)))/sum(taxaIn[,sampleNumber]*as.double(lobo_v)) #raw value
    LOBO20 <- (6.333*LOBO)-5.333
    lobo.results[sampleNumber, ] <- c(LOBO, LOBO20,LOBOtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------LOBO INDEX: END--------############
  #PRECISION
  resultsPath <- resultLoad[[4]]
  precisionmatrix <- read.csv(paste(resultsPath,"\\Precision.csv", sep=""))
  precisionmatrix <- cbind(precisionmatrix, lobo.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="lobo.results$Precision"] <- "LOBO"
  write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$lobo_s > 0)]
  inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv", sep=""))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "LOBO")
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

  rownames(lobo.results) <- resultLoad[[3]]
  return(lobo.results)
}
