#' Calculates the EPID index (EPID)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Dell'Uomo, A. (1996). Assessment of water quality of an Apennine river as a pilot study for diatom-based monitoring of Italian watercourses. Use of algae for monitoring rivers, 65-72.
#' }
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. <https://doi.org/10.1016/j.ecolind.2019.105951>
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
#' epidResults <- diat_epid(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_epid

###### ---------- FUNCTION FOR EPID INDEX (Dell'Uomo) ---------- ########
### INPUT: resultLoad Data. Data needs to be in RA for this index, so if it isn't, the function converts it
### OUTPUTS: dataframe with EPID index per sample
diat_epid <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]] #input data
  #Loads the species list specific for this index
  #epidDB <- read.csv("../Indices/epid.csv") #uses the external csv file
  epidDB <- diathor::epid

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)
  #if acronyms exist, use them, its more precise

  # #exact matches species in input data to acronym from index
  # taxaIn$epid_v <- epidDB$epid_v[match(taxaIn$acronym, trimws(epidDB$acronym))]
  # taxaIn$epid_s <- epidDB$epid_s[match(taxaIn$acronym, trimws(epidDB$acronym))]
  #
  #   #the ones still not found (NA), try against fullspecies
  taxaIn$epid_v <- NA
  taxaIn$epid_s <- NA
    for (i in 1:nrow(taxaIn)) {
      if (is.na(taxaIn$epid_s[i]) | is.na(taxaIn$epid_v[i])){
        taxaIn$epid_v[i] <- epidDB$epid_v[match(trimws(rownames(taxaIn[i,])), trimws(epidDB$fullspecies))]
        taxaIn$epid_s[i] <- epidDB$epid_s[match(trimws(rownames(taxaIn[i,])), trimws(epidDB$fullspecies))]
      }
    }

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------EPID INDEX START --------#############
  print("Calculating EPID index")
  #creates results dataframe
  epid.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(epid.results) <- c("EPID", "EPID20", "Precision")
  #finds the column
  epid_s <- (taxaIn[,"epid_s"])
  epid_v <- (taxaIn[,"epid_v"])
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    EPIDtaxaused <- (length(which(epid_s * taxaIn[,sampleNumber] > 0))*100 / length(epid_s))
    #remove the NA
    epid_s[is.na(epid_s)] = 0
    epid_v[is.na(epid_v)] = 0
    EPID <- sum((taxaIn[,sampleNumber]*as.double(epid_s)*as.double(epid_v)))/sum(taxaIn[,sampleNumber]*as.double(epid_v)) #raw value
    EPID20 <- (-4.75*EPID)+20
    epid.results[sampleNumber, ] <- c(EPID, EPID20,EPIDtaxaused)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------EPID INDEX: END--------############
  #PRECISION
  resultsPath <- resultLoad[[4]]
  #precisionmatrix <- read.csv(paste(resultsPath,"\\Precision.csv", sep=""))
  precisionmatrix <- read.csv(file.path(resultsPath, "Precision.csv"))
  precisionmatrix <- cbind(precisionmatrix, epid.results$Precision)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="epid.results$Precision"] <- "EPID"
  #write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))
  write.csv(precisionmatrix, file.path(resultsPath, "Precision.csv"))
  #END PRECISION



  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$epid_s > 0)]
  #inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv", sep=""))
  inclusionmatrix <- read.csv(file.path(resultsPath, "Taxa included.csv"))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "EPID")
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

  rownames(epid.results) <- resultLoad[[3]]
  return(epid.results)
}

