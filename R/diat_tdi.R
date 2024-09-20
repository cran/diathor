#' Calculates the Trophic (TDI) index
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
#' @description
#' The input for all of these functions is the resulting dataframe (resultLoad) obtained from the diat_loadData() function
#' A CSV or dataframe cannot be used directly with these functions, they have to be loaded first with the diat_loadData() function
#' so the acronyms and species' names are recognized
#' References for the index:
#' \itemize{
#' \item Kelly, M. G., & Whitton, B. A. (1995). The trophic diatom index: a new index for monitoring eutrophication in rivers. Journal of Applied Phycology, 7(4), 433-444.
#' }
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
#' tdiResults <- diat_tdi(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_tdi

###### ---------- FUNCTION FOR TDI INDEX  (Trophic Index - Kelly)  ---------- ########
#### IN THIS SECTION WE CALCULATE TDI INDEX (Trophic Index - Kelly
### INPUT: resultLoad Data cannot be in Relative Abuncance
### OUTPUTS: dataframe with TDI index per sample
diat_tdi <- function(resultLoad, maxDistTaxa = 2){


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
  tdiDB <- diathor::tdi

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  # #the ones still not found (NA), try against fullspecies
  taxaIn$tdi_v <- NA
  taxaIn$tdi_s <- NA


  print("Calculating TDI index")
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$tdi_s[i]) | is.na(taxaIn$tdi_v[i])){
      # New in v0.0.8
      # Uses the stringdist package to find species by names heuristically, with a maximum distance = maxDistTaxa
      # if multiple are found, uses majority consensus to select the correct index value
      # 1) find the species by heuristic search
      spname <- trimws(tolower(rownames(taxaIn[i,])))

      #match the species name to the DB
      species_found <- tdiDB[stringdist::ain(trimws(tolower(tdiDB$fullspecies)),spname, maxDist=maxDistTaxa, matchNA = FALSE),]
      # 2) if found, build majority consensus for sensitivity values
      if (nrow(species_found) > 0){
        vvalue <- as.numeric(names(which.max(table(species_found$tdi_v))))
        svalue <- as.numeric(names(which.max(table(species_found$tdi_s))))
      } else if (nrow(species_found) > 1){
        species_found <- species_found[match(spname, trimws(tolower(species_found$fullspecies)), nomatch=1),]
        vvalue <- as.numeric(names(which.max(table(species_found$tdi_v))))
        svalue <- as.numeric(names(which.max(table(species_found$tdi_s))))
      } else if (nrow(species_found) == 0){
        #species not found, try tautonomy in variety (if species)
        spsplit <- strsplit(spname, " ") #split the name
        #if has epiteth, build new name
        if (length(spsplit[[1]])>1){
          #create vectors with possible epiteths
          newspname <- paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "var.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ") #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "fo.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "f.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "subsp.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "spp.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "ssp.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name
          newspname <- c(newspname, paste(spsplit[[1]][[1]], spsplit[[1]][[2]], "var.", spsplit[[1]][[2]], "fo.", spsplit[[1]][[length(spsplit[[1]])]], sep = " ")) #create new sp name

          #search again against all possible epiteths
          species_found <- tdiDB[stringdist::ain(trimws(tolower(tdiDB$fullspecies)),newspname, maxDist=maxDistTaxa, matchNA = FALSE),]
          if (nrow(species_found) > 0){
            #found with tautonomy
            vvalue <- as.numeric(names(which.max(table(species_found$tdi_v[1]))))
            svalue <- as.numeric(names(which.max(table(species_found$tdi_s[1]))))
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
    }
      #records the final consensus value
      taxaIn$tdi_v[i] <- vvalue
      taxaIn$tdi_s[i] <- svalue
    }



  #######--------TDI INDEX START --------#############

  #creates results dataframe
  tdi.results <- data.frame(matrix(ncol = 3, nrow = (lastcol-1)))
  colnames(tdi.results) <- c("TDI20", "TDI100", "num_taxa")
  #finds the column
  tdi_s <- (taxaIn[,"tdi_s"])
  tdi_v <- (taxaIn[,"tdi_v"])

  # Prints the number of taxa recognized for this index, regardless of their abundance
  # It is therefore the same for all samples

  number_recognized_taxa <- round((100 - (sum(is.na(taxaIn$tdi_s)) / nrow(taxaIn))*100),1)
  print(paste("Taxa recognized to be used in TDI index: ", number_recognized_taxa, "%"))




  #PROGRESS BAR
  if (interactive()) {
    pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  }

  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #how many taxa will be used to calculate?
    #Revised v0.0.8
    num_taxa <- length(which(tdi_s * taxaIn[,sampleNumber] > 0))
    #remove the NA
    tdi_s[is.na(tdi_s)] = 0
    tdi_v[is.na(tdi_v)] = 0
    TDI <- sum((taxaIn[,sampleNumber]*as.double(tdi_s)*as.double(tdi_v)))/sum(taxaIn[,sampleNumber]*as.double(tdi_v)) #raw value
    TDI20 <- (-4.75*TDI)+24.75
    TDI100 <- (TDI*25)-25
    tdi.results[sampleNumber, ] <- c(TDI20, TDI100,num_taxa)
    #update progressbar
    if (interactive()) {
      setTxtProgressBar(pb, sampleNumber)
    }
  }
  #close progressbar
  if (interactive()) {
    close(pb)
  }
  #######--------TDI INDEX: END--------############

  #PRECISION RECORDING
  resultsPath <- resultLoad[[4]]
  #reads the csv file
  precisionmatrix <- read.csv(file.path(resultsPath, "num_taxa.csv"))
  #joins with the precision column
  precisionmatrix <- cbind(precisionmatrix, tdi.results$num_taxa)
  precisionmatrix <- precisionmatrix[-(1:which(colnames(precisionmatrix)=="Sample")-1)]
  names(precisionmatrix)[names(precisionmatrix)=="tdi.results$num_taxa"] <- "TDI"
  write.csv(precisionmatrix, file.path(resultsPath, "num_taxa.csv"))
  #END PRECISION

  #TAXA INCLUSION
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(taxaIn$tdi_s > 0)]
  #inclusionmatrix <- read.csv(paste(resultsPath,"\\Taxa included.csv", sep=""))
  inclusionmatrix <- read.csv(file.path(resultsPath, "Taxa included.csv"))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "TDI")
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
  #write.csv(exclusionmatrix, paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  write.csv(exclusionmatrix, file.path(resultsPath,"Taxa excluded.csv"))
  #END EXCLUDED TAXA

  rownames(tdi.results) <- resultLoad[[3]]
  return(tdi.results)
}

