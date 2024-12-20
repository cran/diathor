#' Calculates the Swedish Phosphorus Diatom Index (PDIse)
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Kahlert, M., Fölster, J., & Tapolczai, K. (2023). No lukewarm diatom communities—the response of freshwater benthic diatoms to phosphorus in streams as basis for a new phosphorus diatom index (PDISE). Environmental Monitoring and Assessment, 195(7), 807.
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
#' pdiseResults <- diat_pdise(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_pdise


###### ---------- FUNCTION FOR PDIse INDEX  (Kahlert et al. 2023)---------- ########
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with PDIse index per sample

diat_pdise <- function(resultLoad, maxDistTaxa = 2){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculation of the PDIse cancelled")
    }
  }

  #pre-test to check if data is abundance data. For the PDIse, data has to be abundance data, it cannot be relative abundance!!!
  taxaIn <- resultLoad[[2]] #supposed to be real abundance
  new_species_vec <- taxaIn$new_species
  numeric_cols <- taxaIn[sapply(taxaIn, is.numeric)]
  if (all(colSums(numeric_cols) < 101)) {
    print(paste("Data for the PDIse has to be absolute abundance data, not relative abundance"))
    stop("Calculation of the PDIse cancelled")
  } else {
    #If all columns add to 101 or less (dont use 100 since sometimes decimals end up adding to a bit more than 100)
    #then convert the data to square root transformed abundance (as specified in the original paper)
    taxaIn <- sqrt(numeric_cols)
    taxaIn$new_species <- new_species_vec
  }

  #Loads the species list specific for this index
  pdiseDB <- diathor::pdise

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  # #exact matches species in input data to acronym from index

  taxaIn$pdise_v <- NA
  taxaIn$pdise_s <- NA

  print("Calculating PDIse index - BETA (tread lightly)")
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$pdise_v[i]) | is.na(taxaIn$pdise_s[i])){
      # New from v0.0.8 onwards
      # Uses the stringdist package to find species by names heuristically, with a maximum distance = maxDistTaxa
      # if multiple are found, uses majority consensus to select the correct index value
      # 1) find the species by heuristic search.
      spname <- trimws(tolower(rownames(taxaIn[i,])))

      species_found <- pdiseDB[stringdist::ain(trimws(tolower(pdiseDB$fullspecies)),spname, maxDist=maxDistTaxa, matchNA = FALSE),]
      # 2) if found, build majority consensus for sensitivity values
      if (nrow(species_found) == 1){
        vvalue <- as.numeric(names(which.max(table(species_found$pdise_v))))
        svalue <- as.numeric(names(which.max(table(species_found$pdise_s))))
        taxaIn$new_species[i] <- species_found$fullspecies[1]
      } else if (nrow(species_found) > 1){
        species_found <- species_found[match(spname, trimws(tolower(species_found$fullspecies)), nomatch=1),]
        vvalue <- as.numeric(names(which.max(table(species_found$pdise_v))))
        svalue <- as.numeric(names(which.max(table(species_found$pdise_s))))
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
          species_found <- pdiseDB[stringdist::ain(trimws(tolower(pdiseDB$fullspecies)),newspname, maxDist=maxDistTaxa, matchNA = FALSE),]
          if (nrow(species_found) > 0){
            #found with tautonomy
            vvalue <- as.numeric(names(which.max(table(species_found$pdise_v[1]))))
            svalue <- as.numeric(names(which.max(table(species_found$pdise_s[1]))))
            taxaIn$new_species[i] <- species_found$fullspecies[1]
          } else {
            #species not found, make everything NA
            vvalue = NA
            svalue = NA
          }
        }  else {
          #species not found, make everything NA
          vvalue = NA
          svalue = NA
        }
      }
      #records the final consensus value
      taxaIn$pdise_v[i] <- vvalue
      taxaIn$pdise_s[i] <- svalue
    }
  }

  #removes NA from taxaInRA
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------PDIse INDEX START --------#############

  #creates results dataframe
  pdise.results <- data.frame(matrix(ncol = 2, nrow = (lastcol-1)))
  colnames(pdise.results) <- c("PDIse", "num_taxa")
  #finds the column
  pdise_s <- (taxaIn[,"pdise_s"])
  pdise_v <- (taxaIn[,"pdise_v"])

  # Prints the number of taxa recognized for this index, regardless of their abundance
  # It is therefore the same for all samples

  number_recognized_taxa <- round((100 - (sum(is.na(taxaIn$pdise_s)) / nrow(taxaIn))*100),1)
  print(paste("Taxa recognized to be used in PDIse index: ", number_recognized_taxa, "%"))


  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    # Added in v0.0.8
    num_taxa <- length(which(pdise_s * taxaIn[,sampleNumber] > 0))
    #remove the NA
    pdise_s[is.na(pdise_s)] = 0
    pdise_v[is.na(pdise_v)] = 0
    PDIse <- sum((taxaIn[,sampleNumber]*as.double(pdise_s)*as.double(pdise_v)))/sum(taxaIn[,sampleNumber]*as.double(pdise_v)) #raw value
    pdise.results[sampleNumber, ] <- c(PDIse, num_taxa)
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)

  #######--------PDIse INDEX: END--------############
  #PRECISION RECORDING
  resultsPath <- resultLoad[[4]]
  #reads the csv file
  precisionmatrix <- read.csv(file.path(resultsPath, "num_taxa.csv"))
  #joins with the precision column
  precisionmatrix <- cbind(precisionmatrix, pdise.results$num_taxa)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="pdise.results$num_taxa"] <- "PDIse"
  write.csv(precisionmatrix, file.path(resultsPath, "num_taxa.csv"))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$pdise_s > 0)]
  inclusionmatrix <- read.csv(file.path(resultsPath, "Taxa included.csv"))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "PDIse")
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
  } else{newinclusionmatrix[is.na(newinclusionmatrix) == FALSE] <- NA}

  inclusionmatrix <- newinclusionmatrix
  colnames(inclusionmatrix) <- colnamesInclusionMatrix
  inclusionmatrix <- inclusionmatrix[-(1:which(colnames(inclusionmatrix)=="Eco.Morpho")-1)]
  write.csv(inclusionmatrix, file.path(resultsPath,"Taxa included.csv"))
  #END TAXA INCLUSION
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
  }else{newexclusionmatrix[is.na(newexclusionmatrix) == FALSE] <- NA}

  exclusionmatrix <- newexclusionmatrix
  colnames(exclusionmatrix) <- colnamesInclusionMatrix
  exclusionmatrix <- exclusionmatrix[-(1:which(colnames(exclusionmatrix)=="Eco.Morpho")-1)]
  write.csv(exclusionmatrix, file.path(resultsPath,"Taxa excluded.csv"))
  #END EXCLUDED TAXA

  rownames(pdise.results) <- resultLoad[[3]]
  return(pdise.results)
}

