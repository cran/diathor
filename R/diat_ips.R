#' Calculates the Specific Polluosensitivity Index (IPS) index
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Coste, M. (1982). Étude des méthodes biologiques d’appréciation quantitative de la qualité des eaux. Rapport Cemagref QE Lyon-AF Bassin Rhône Méditerranée Corse.
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
#' ipsResults <- diat_ips(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_ips


###### ---------- FUNCTION FOR IPS INDEX   ---------- ########
#### IN THIS SECTION WE CALCULATE IPS INDEX (REF)
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with IPS index per sample
diat_ips <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #Loads the species list specific for this index
  #idpDB <- read.csv("../Indices/ips.csv") #uses the external csv file
  ipsDB <- diathor::ips

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  #if acronyms exist, use them, its more precise
  #if there is an acronym column, it removes it and stores it for later
  # #exact matches species in input data to acronym from index
  # taxaIn$ips_s <- ipsDB$ips_s[match(trimws(taxaIn$acronym), trimws(ipsDB$acronym))]
  # taxaIn$ips_v <- ipsDB$ips_v[match(trimws(taxaIn$acronym), trimws(ipsDB$acronym))]


  # #the ones still not found (NA), try against fullspecies
  taxaIn$ips_v <- NA
  taxaIn$ips_s <- NA
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$ips_s[i]) | is.na(taxaIn$ips_v[i])){
      taxaIn$ips_v[i] <- ipsDB$ips_v[match(trimws(rownames(taxaIn[i,])), trimws(ipsDB$fullspecies))]
      taxaIn$ips_s[i] <- ipsDB$ips_s[match(trimws(rownames(taxaIn[i,])), trimws(ipsDB$fullspecies))]
    }
  }


  #######--------IPS INDEX START (indice poluto sensible)--------#############
  print("Calculating IPS index")
  #creates results dataframe
  ips.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(ips.results) <- c("IPS", "IPS20", "Precision")
  #finds the column
  ips_s <- (taxaIn[,"ips_s"])
  ips_v <- (taxaIn[,"ips_v"])

  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    IPStaxaused <- (length(which(ips_s * taxaIn[,sampleNumber] > 0))*100 / length(ips_s))
    #remove the NA
    ips_s[is.na(ips_s)] = 0
    ips_v[is.na(ips_v)] = 0
    IPS <- sum((taxaIn[,sampleNumber]*as.double(ips_s)*as.double(ips_v)))/sum(taxaIn[,sampleNumber]*as.double(ips_v)) #raw value
    IPS20 <- (IPS*4.75)-3.75 #STANDARDIZED VALUE TO 20
    ips.results[sampleNumber, ] <- c(IPS, IPS20, IPStaxaused)

    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------IPS INDEX: END--------############

  #PRECISION
  resultsPath <- resultLoad[[4]]
  precisionmatrix <- read.csv(paste(resultsPath,"\\Precision.csv", sep=""))
  precisionmatrix <- cbind(precisionmatrix, ips.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="ips.results$Precision"] <- "IPS"
  write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$ips_s > 0)]
  inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv", sep=""))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "IPS")
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
  #END INCLUDED TAXA
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

  rownames(ips.results) <- resultLoad[[3]]
  return(ips.results)
}

