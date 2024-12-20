#' Calculates the Specific Polluosensitivity Index (IPS) index
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
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
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi:10.1016/j.ecolind.2019.105951
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
diat_ips <- function(resultLoad, maxDistTaxa = 2){
  #DEBUG: resulting time check
  # start.time <- Sys.time()

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }
  #Gets the resultLoad from the diatLoad function
  taxaIn <- resultLoad[[2]]

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #Loads the species list specific for this index
  ipsDB <- diathor::ips

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  # Attempts to match full species names
  taxaIn$ips_v <- NA
  taxaIn$ips_s <- NA
  print("Calculating IPS index")
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$ips_s[i]) | is.na(taxaIn$ips_v[i])){
      # New in v0.0.8
      # Uses the stringdist package to find species by names heuristically, with a maximum distance = maxDistTaxa
      # if multiple are found, uses majority consensus to select the correct index value
      # 1) find the species by heuristic search
      spname <- trimws(tolower(rownames(taxaIn[i,])))

      species_found <- ipsDB[stringdist::ain(trimws(tolower(ipsDB$fullspecies)),spname, maxDist=maxDistTaxa, matchNA = FALSE),]
      # 2) if found, build majority consensus for sensitivity values
      if (nrow(species_found) == 1){
        vvalue <- as.numeric(names(which.max(table(species_found$ips_v))))
        svalue <- as.numeric(names(which.max(table(species_found$ips_s))))
        taxaIn$new_species[i] <- species_found$fullspecies[1]
      } else if (nrow(species_found) > 1){
        species_found <- species_found[match(spname, trimws(tolower(species_found$fullspecies)), nomatch=1),]
        vvalue <- as.numeric(names(which.max(table(species_found$ips_v))))
        svalue <- as.numeric(names(which.max(table(species_found$ips_s))))
      } else if (nrow(species_found) == 0){
        #species not found, try tautonomy in variety
        spsplit <- strsplit(spname, " ") #split the name
        #if has epiteth
        if (length(spsplit[[1]])>1){
          #create vectors with possible epiteths
          newspname <- paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "var.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ") #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "fo.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "subsp.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "spp.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "ssp.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "var.", spsplit[[1]][[2]], "fo.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name

          #search again against all possible epiteths
          species_found <- ipsDB[stringdist::ain(trimws(tolower(ipsDB$fullspecies)),newspname, maxDist=maxDistTaxa, matchNA = FALSE),]
          if (nrow(species_found) > 0){
            #found with tautonomy
            vvalue <- as.numeric(names(which.max(table(species_found$ips_v[1]))))
            svalue <- as.numeric(names(which.max(table(species_found$ips_s[1]))))
            taxaIn$new_species[i] <- species_found$fullspecies[1]
          } else {
              #species not found, make everything NA
              vvalue = NA
              svalue = NA
          }
        }  else {
          # length(spsplit[[1]]) =<1
          #species not found, make everything NA
          vvalue = NA
          svalue = NA
        }
      }
      #records the final consensus value
      taxaIn$ips_v[i] <- vvalue
      taxaIn$ips_s[i] <- svalue
    }
  }




  #######--------IPS INDEX START (indice poluto sensible)--------#############

  #creates results dataframe
  ips.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(ips.results) <- c("IPS", "IPS20", "num_taxa")
  #finds the column
  ips_s <- (taxaIn$ips_s)
  ips_v <- (taxaIn$ips_v)

  # Prints the number of taxa recognized for this index, regardless of their abundance
  # It is therefore the same for all samples

  number_recognized_taxa <- round((100 - (sum(is.na(taxaIn$ips_s)) / nrow(taxaIn))*100),1)
  print(paste("Taxa recognized to be used in IPS index: ", number_recognized_taxa, "%"))


  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate with abundance > 0? (precision)
    #Revised v0.0.8
    num_taxa <- length(which(ips_s * taxaIn[,sampleNumber] > 0))

    #remove the NA
    ips_s[is.na(ips_s)] = 0
    ips_v[is.na(ips_v)] = 0
    IPS <- sum((taxaIn[,sampleNumber]*as.double(ips_s)*as.double(ips_v)))/sum(taxaIn[,sampleNumber]*as.double(ips_v)) #raw value
    IPS20 <- (IPS*4.75)-3.75 #STANDARDIZED VALUE TO 20
    #Record values
    ips.results[sampleNumber, ] <- c(IPS, IPS20, num_taxa)

    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  #######--------IPS INDEX: END--------############

  #PRECISION RECORDING
  resultsPath <- resultLoad[[4]]
  #reads the csv file
  precisionmatrix <- read.csv(file.path(resultsPath, "num_taxa.csv"))
  #joins with the precision column
  precisionmatrix <- cbind(precisionmatrix, ips.results$num_taxa)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="ips.results$num_taxa"] <- "IPS"
  #writes the precision file
  write.csv(precisionmatrix, file.path(resultsPath, "num_taxa.csv"))
  #END PRECISION RECORDING

  #TAXA INCLUSION
  taxaIncluded <- taxaIn$species[which(taxaIn$ips_s > 0)]
  inclusionmatrix <- read.csv(file.path(resultsPath, "Taxa included.csv"))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "IPS")
  #creates a new data matrix to append the existing Taxa Included file
  newinclusionmatrix <- as.data.frame(matrix(nrow=max(length(taxaIncluded), nrow(inclusionmatrix)), ncol=ncol(inclusionmatrix)+1))
  for (i in 1:ncol(inclusionmatrix)){
    newinclusionmatrix[1:nrow(inclusionmatrix),i] <- as.character(inclusionmatrix[1:nrow(inclusionmatrix),i])
  }
  #check that taxaIncluded is at least 1
  if (length(taxaIncluded) > 0) {
    if (nrow(newinclusionmatrix) > length(taxaIncluded)){
      newinclusionmatrix[1:length(taxaIncluded), ncol(newinclusionmatrix)] <- taxaIncluded
    } else {
      newinclusionmatrix[1:nrow(newinclusionmatrix), ncol(newinclusionmatrix)] <- taxaIncluded
    }
  }else{newinclusionmatrix[is.na(newinclusionmatrix) == FALSE] <- NA}

  inclusionmatrix <- newinclusionmatrix
  colnames(inclusionmatrix) <- colnamesInclusionMatrix
  inclusionmatrix <- inclusionmatrix[-(1:which(colnames(inclusionmatrix)=="Eco.Morpho")-1)]
  write.csv(inclusionmatrix, file.path(resultsPath,"Taxa included.csv"))
  #END INCLUDED TAXA

  #EXCLUDED TAXA
  taxaExcluded <- taxaIn[!('%in%'(taxaIn$species,taxaIncluded)),"species"]
  exclusionmatrix <- read.csv(file.path(resultsPath, "Taxa excluded.csv"))
  #creates a new data matrix to append the existing Taxa Included file
  newexclusionmatrix <- as.data.frame(matrix(nrow=max(length(taxaExcluded), nrow(exclusionmatrix)), ncol=ncol(exclusionmatrix)+1))
  for (i in 1:ncol(exclusionmatrix)){
    newexclusionmatrix[1:nrow(exclusionmatrix),i] <- as.character(exclusionmatrix[1:nrow(exclusionmatrix),i])
  }
  #check that taxaExcluded is at least 1
  if (length(taxaExcluded) > 0) {
    if (nrow(newexclusionmatrix) > length(taxaExcluded)){
      newexclusionmatrix[1:length(taxaExcluded), ncol(newexclusionmatrix)] <- taxaExcluded
    } else {
      newexclusionmatrix[1:nrow(newexclusionmatrix), ncol(newexclusionmatrix)] <- taxaExcluded
    }
  } else{newexclusionmatrix[is.na(newexclusionmatrix) == FALSE] <- NA}
  exclusionmatrix <- newexclusionmatrix
  colnames(exclusionmatrix) <- colnamesInclusionMatrix
  exclusionmatrix <- exclusionmatrix[-(1:which(colnames(exclusionmatrix)=="Eco.Morpho")-1)]
  #write.csv(exclusionmatrix, paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  write.csv(exclusionmatrix, file.path(resultsPath,"Taxa excluded.csv"))
  #END EXCLUDED TAXA

  rownames(ips.results) <- resultLoad[[3]]

  #DEBUG: printint calculation timer
  # end.time <- Sys.time()
  # time.taken <- end.time - start.time
  # print(paste("time", time.taken))

  return(ips.results)
}

