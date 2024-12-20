#' Calculate the combined classification of ecological guilds and size classes for diatoms
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for these functions is the resulting dataframe obtained from the diat_loadData() function, to calculate the ecological guilds for the diatoms
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi:10.1016/j.ecolind.2019.105951
#' }
#' Classification is obtained from:
#' \itemize{
#' \item B-Béres, V., Török, P., Kókai, Z., Lukács, Á., Enikő, T., Tóthmérész, B., & Bácsi, I. (2017). Ecological background of diatom functional groups: Comparability of classification systems. Ecological Indicators, 82, 183-188.
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
#' guildsResults <- diat_cemfgs_rb(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_cemfgs_rb
#'

###### ---------- FUNCTION FOR ECOLOGICAL GUILDS & SIZE CLASSES   ---------- ########
### INPUT: taxaInRA created in loadData()
### OUTPUTS: dataframe with ecological guilds  per sample

diat_cemfgs_rb <- function(resultLoad){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[1]] #1 = Relative abundance, 2=abundance

  #Loads the species list specific for this index
  cemfgs_rbDB <- diathor::cemfgs_rb

  #creates a species column with the rownames to fit in the script
  taxaIn$species <- row.names(taxaIn)

  # try finding fullspecies
  taxaIn$cemfgs_rb <- NA
  for (i in 1:nrow(taxaIn)) {
    if (is.na(taxaIn$cemfgs_rb[i])){
      taxaIn$cemfgs_rb[i] <- cemfgs_rbDB$cemfgs_rb[match(trimws(rownames(taxaIn[i,])), trimws(cemfgs_rbDB$fullspecies))]
    }
  }

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #######--------CEMFGS_RG INDEX START --------#############
  print("Calculating CEMFGS_RG")
  #creates results dataframe
  cemfgs_rb.results <- data.frame(matrix(ncol = 2, nrow = (lastcol-1)))
  colnames(cemfgs_rb.results) <- c("CEMFGS_RG", "Precision")
  #finds the column
  cemfgs_rb <- (taxaIn[,"cemfgs_rb"])

  #PROGRESS BAR
  # pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  HS1 <- HS2 <- HS3 <- HS4 <- HS5 <- NULL
  LS1 <- LS2 <- LS3 <- LS4 <- LS5 <- NULL
  MS1 <- MS2 <- MS3 <- MS4 <- MS5 <- NULL
  PS1 <- PS2 <- PS3 <- PS4 <- PS5 <- NULL
  CEMFGS_class_Indet <- CEMFGS_Taxa_used <- CEMFGS_RB_Indet <- NULL


  #OBSOLETE CODE FROM V0.1.4
  # data.table::setDT(taxaIn)
  # cemfgs_rb.results <- suppressWarnings(data.table(
  #   HS1 <- unlist(taxaIn[which(cemfgs_rb == "HS1"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                            1)]),
  #   HS2 <- unlist(taxaIn[which(cemfgs_rb == "HS2"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   HS3 <- unlist(taxaIn[which(cemfgs_rb == "HS3"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   HS4 <- unlist(taxaIn[which(cemfgs_rb == "HS4"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   HS5 <- unlist(taxaIn[which(cemfgs_rb == "HS5"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   LS1 <- unlist(taxaIn[which(cemfgs_rb == "LS1"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   LS2 <- unlist(taxaIn[which(cemfgs_rb == "LS2"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   LS3 <- unlist(taxaIn[which(cemfgs_rb == "LS3"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   LS4 <- unlist(taxaIn[which(cemfgs_rb == "LS4"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   LS5 <- unlist(taxaIn[which(cemfgs_rb == "LS5"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   MS1 <- unlist(taxaIn[which(cemfgs_rb == "MS1"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   MS2 <- unlist(taxaIn[which(cemfgs_rb == "MS2"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   MS3 <- unlist(taxaIn[which(cemfgs_rb == "MS3"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   MS4 <- unlist(taxaIn[which(cemfgs_rb == "MS4"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   MS5 <- unlist(taxaIn[which(cemfgs_rb == "MS5"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   PS1 <- unlist(taxaIn[which(cemfgs_rb == "PS1"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   PS2 <- unlist(taxaIn[which(cemfgs_rb == "PS2"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   PS3 <- unlist(taxaIn[which(cemfgs_rb == "PS3"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   PS4 <- unlist(taxaIn[which(cemfgs_rb == "PS4"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)]),
  #   PS5 <- unlist(taxaIn[which(cemfgs_rb == "PS5"), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
  #                                                                                                 1)])
  # ))
  #
  #

  ### NEW CODE FROM V0.1.5 ONWARDS
  data.table::setDT(taxaIn)

  # Define the groups
  groups <- c("HS1", "HS2", "HS3", "HS4", "HS5",
              "LS1", "LS2", "LS3", "LS4", "LS5",
              "MS1", "MS2", "MS3", "MS4", "MS5",
              "PS1", "PS2", "PS3", "PS4", "PS5")

  # Initialize an empty data.table
  cemfgs_rb.results <- data.table(matrix(0, nrow = (lastcol - 1), ncol = length(groups)))

  # Add column names
  setnames(cemfgs_rb.results, groups)

  # Calculate sums and add to results
  for (group in groups) {
    group_data <- taxaIn[cemfgs_rb == group,
                         lapply(.SD, sum, na.rm = TRUE),
                         .SDcols = 1:(lastcol - 1)]

    if (nrow(group_data) > 0) {
      cemfgs_rb.results[[group]] <- unlist(group_data)
    }
  }


  #### END NEW VERSION

    #replace NAs for 0
  cemfgs_rb.results[is.na(cemfgs_rb.results)] = 0
  cemfgs_rb.results[, `CEMFGS_RB_Indet` := round(100 - ( HS1 + HS2 + HS3 + HS4 + HS5 + LS1 + LS2 + LS3 + LS4 + LS5 + MS1 + MS2 + MS3 + MS4 + MS5 + PS1 + PS2 + PS3 + PS4 + PS5), 1)]


  #taxa used for each guild
  cemfgs_rb_binary <- as.numeric(!is.na(taxaIn$cemfgs_rb))

  for (sampleNumber in 1:(lastcol - 1)) {
    #CEMFGS_RGtaxaused <- length(which(cemfgs_rb_binary * taxaIn[,..sampleNumber] > 0))
    CEMFGS_RGtaxaused <- length(which(cemfgs_rb_binary * taxaIn[,sampleNumber, with = F] > 0))

    #CEMFGS_RGtaxaused <- length(which(cemfgs_rb_binary * taxaIn[,..sampleNumber] > 0))*100 / length(cemfgs_rb_binary)
    #The .. before sampleNumber are the new way to reference the variable in the data.table package
    cemfgs_rb.results$CEMFGS_RB_Taxa_used [sampleNumber] <- CEMFGS_RGtaxaused

  }

  cemfgs_rb.results2 <- matrix(unlist(cemfgs_rb.results), ncol = lastcol - 1, byrow = TRUE) %>% as.data.frame()
  cemfgs_rb.results2 <- t(cemfgs_rb.results2)
  rownames(cemfgs_rb.results2) <- resultLoad[[3]]
  colnames(cemfgs_rb.results2) <- c("HS1" , "HS2" , "HS3" , "HS4" , "HS5" , "LS1" , "LS2" , "LS3" , "LS4" , "LS5" , "MS1" , "MS2" , "MS3" , "MS4" , "MS5" , "PS1" , "PS2" , "PS3" , "PS4" , "PS5" , "CEMFGS class Indet" , "CEMFGS Taxa used")

  cemfgs_rb.results <- as.data.frame(cemfgs_rb.results2) #need to convert it to dataframe explicitly to plot

  #TAXA INCLUSION
  resultsPath <- resultLoad[[4]]
  #taxa with acronyms
  taxaIncluded <- taxaIn$species[which(!is.na(taxaIn$cemfgs_rb))]
  inclusionmatrix <- read.csv(file.path(resultsPath, "Taxa included.csv"))
  colnamesInclusionMatrix <- c(colnames(inclusionmatrix), "CEMFGS_RB")
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
    if (nrow(newexclusionmatrix) > nrow(taxaExcluded)){
      newexclusionmatrix[1:nrow(taxaExcluded), ncol(newexclusionmatrix)] <- taxaExcluded
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






  return(cemfgs_rb.results)

}

