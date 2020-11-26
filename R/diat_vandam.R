#' Calculates ecological information for diatoms based on the Van Dam classification
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @param vandamReports Boolean. If set to 'TRUE' the detailed reports for the Van Dam classifications will be reported in the Output. Default = TRUE
#' @description
#' The input for these functions is the resulting dataframe obtained from the diat_loadData() function, to calculate ecological information for diatoms based on the Van Dam classification
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
#' }
#' Van Dam classification is obtained form:
#' \itemize{
#' \item Van Dam, H., Mertens, A., & Sinkeldam, J. (1994). A coded checklist and ecological indicator values of freshwater diatoms from the Netherlands. Netherland Journal of Aquatic Ecology, 28(1), 117-133.
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
#' vandamResults <- diat_vandam(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_vandam
#'

###### ---------- FUNCTION FOR VANDAM ECOLOGICAL CHARACTERIZATION: DONE   ---------- ########
#### IN THIS SECTION WE CALCULATE ECOLOGICAL PREFERENCES ACCORDING TO VAN DAM (REF)
### INPUT: resultLoad created in loadData()
### OUTPUTS: dataframe with VanDam's ecological values
### PARAMETERS: vandamReports (Boolean), if a detailed report of the taxa used for each sub index is exported or not
#vdams = salinity
#vdamnh = N-Heterotrophie
#vdamo2 = Oxygen
#vdamsap = saprobity
#vdamtrop = trophic
#vdamaero = aerophilie

diat_vandam <- function(resultLoad, vandamReports=TRUE){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaInEco <- resultLoad[[5]]
  #checks thata taxaInEco (taxaInEco from diat_loadData) has at least recognized some species
  if (nrow(taxaInEco)==0){
    print("No species were recognized for VanDam calculations")
    print("VanDam data will not be available")
    vandam.results <- NULL
    return(vandam.results)
  }

  #gets the column named "species", everything before that is a sample
  lastcol <- which(colnames(taxaInEco)=="species")

  #Convert taxaIn sample data to Relative Abundance data
  taxaInRA <- taxaInEco
  for (i in 1:nrow(taxaInEco)){
    for (j in 1:(lastcol-1)){
      if (is.na(taxaInEco[i,j])){
        taxaInRA[i,j] <- 0
      } else {
        taxaInRA[i,j] <- (taxaInEco[i,j]*100)/sum(taxaInEco[,j])
      }
    }
  }

  #removes NA from taxaInRA
  taxaInRA[is.na(taxaInRA)] <- 0

  resultsPath <- resultLoad[[4]]
  # CHECKS IF RESULTSPATH EXISTS, OR ASKS FOR IT, FOR THE DETAILED REPORTS
  if (vandamReports==TRUE & is.na(resultsPath)==TRUE){
    print("Select Results folder for the detailed VanDam reports")
    resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  }

  #EXPORT REPORTS FOR EACH SAMPLE OF TAXA USED
  if (vandamReports==TRUE & is.na(resultsPath)==FALSE){
    print("Exporting detailed reports for VanDam ecological preferences")
    print(resultsPath)
    vandamtaxafile = paste("VanDam Taxa used.txt", sep="")
    write("TAXA USED FOR EACH ECOLOGICAL VARIABLE USING VANDAM's CLASSIFICATION", paste(resultsPath, "\\", vandamtaxafile, sep=""))
    write("These taxa were included because: ",paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)
    write("a) they had a reliable classification in VanDam's classification system",paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)
    write("b) they had a relative abundance in the sample > 0",paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)
  }

  ##### VANDAM SALINITY
  #creates results dataframe
  vdams_labels <- c("VD Salinity 1", "VD Salinity 2", "VD Salinity 3", "VD Salinity 4", "VD Salinity Indet", "VD Salinity Taxa used")
  vdamSalinity <- data.frame(matrix(ncol = 6, nrow = (lastcol-1)))
  colnames(vdamSalinity) <- vdams_labels

  #finds the column
  vdams_v <- (taxaInRA[,"vdams"])
  #remove the NA
  vdams_v[is.na(vdams_v)] = 0
  print("Calculating Van Dam salinity")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdams_cat1 <- sum(taxaInRA[which(vdams_v == 1),sampleNumber])
    vdams_cat2 <- sum(taxaInRA[which(vdams_v == 2),sampleNumber])
    vdams_cat3 <- sum(taxaInRA[which(vdams_v == 3),sampleNumber])
    vdams_cat4 <- sum(taxaInRA[which(vdams_v == 4),sampleNumber])
    vdams_indet <- 100 - sum(vdams_cat1, vdams_cat2, vdams_cat3, vdams_cat4) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    Vdamstaxaused <- length(which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamstaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    vdams_values <- c(vdams_cat1, vdams_cat2, vdams_cat3, vdams_cat4, vdams_indet, Vdamstaxaused)
    vdamSalinity[sampleNumber, ] <- vdams_values
    if (vandamReports==TRUE & exists("resultsPath")) {write(paste("SALINITY- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
    if (vandamReports==TRUE & exists("resultsPath")) {write(Vdamstaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
  }

  ##### END VANDAM SALINITY

  ##### VANDAM N-HETEROTROPHIE
  #creates results dataframe
  vdamnh_labels <- c("VD N-Het 1", "VD N-Het 2", "VD N-Het 3", "VD N-Het 4", "VD N-Het Indet", "VD N-Het Taxa used")
  vdamNHeterotrophy <- data.frame(matrix(ncol = 6, nrow = (lastcol-1)))
  colnames(vdamNHeterotrophy) <- vdamnh_labels
  #finds the column
  vdamnh_v <- (taxaInRA[,"vdamnh"])
  #remove the NA
  vdamnh_v[is.na(vdamnh_v)] = 0
  print("Calculating Van Dam N-heterotrophy")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamnh_cat1 <- sum(taxaInRA[which(vdamnh_v == 1),sampleNumber])
    vdamnh_cat2 <- sum(taxaInRA[which(vdamnh_v == 2),sampleNumber])
    vdamnh_cat3 <- sum(taxaInRA[which(vdamnh_v == 3),sampleNumber])
    vdamnh_cat4 <- sum(taxaInRA[which(vdamnh_v == 4),sampleNumber])
    vdamnh_indet <- 100 - sum(vdamnh_cat1, vdamnh_cat2, vdamnh_cat3, vdamnh_cat4) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    Vdamnhtaxaused <- length(which(vdamnh_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamnhtaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    vdamnh_values <- c(vdamnh_cat1, vdamnh_cat2, vdamnh_cat3, vdamnh_cat4, vdamnh_indet, Vdamnhtaxaused)
    vdamNHeterotrophy[sampleNumber, ] <- vdamnh_values
    if (vandamReports==TRUE & exists("resultsPath")) {write(paste("N-HETEROTROPHY- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
    if (vandamReports==TRUE & exists("resultsPath")) {write(Vdamnhtaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
  }
  ##### END VANDAM N-HETEROTROPHIE

  ##### VANDAM OXYGEN REQUIREMENTS
  #creates results dataframe
  vdamo2_labels <- c("VD Oxygen 1", "VD Oxygen 2", "VD Oxygen 3", "VD Oxygen 4","VD Oxygen 5", "VD Oxygen Indet", "VD Oxygen Taxa used")
  vdamOxygen <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(vdamOxygen) <- vdamo2_labels
  #finds the column
  vdamo2_v <- (taxaInRA[,"vdamo2"])
  #remove the NA
  vdamo2_v[is.na(vdamo2_v)] = 0
  print("Calculating Van Dam oxygen requirements")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamo2_cat1 <- sum(taxaInRA[which(vdamo2_v == 1),sampleNumber])
    vdamo2_cat2 <- sum(taxaInRA[which(vdamo2_v == 2),sampleNumber])
    vdamo2_cat3 <- sum(taxaInRA[which(vdamo2_v == 3),sampleNumber])
    vdamo2_cat4 <- sum(taxaInRA[which(vdamo2_v == 4),sampleNumber])
    vdamo2_cat5 <- sum(taxaInRA[which(vdamo2_v == 5),sampleNumber])
    vdamo2_indet <- 100 - sum(vdamo2_cat1, vdamo2_cat2, vdamo2_cat3, vdamo2_cat4, vdamo2_cat5) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    Vdamo2taxaused <- length(which(vdamo2_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamo2taxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    vdamo2_values <- c(vdamo2_cat1, vdamo2_cat2, vdamo2_cat3, vdamo2_cat4, vdamo2_cat5, vdamo2_indet, Vdamo2taxaused)
    vdamOxygen[sampleNumber, ] <- vdamo2_values
    if (vandamReports==TRUE & exists("resultsPath")) {write(paste("OXYGEN REQUIREMENTS- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
    if (vandamReports==TRUE & exists("resultsPath")) {write(Vdamo2taxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
  }
  ##### END VANDAM OXYGEN REQUIREMENTS

  ##### VANDAM SAPROBITY
  #creates results dataframe
  vdamsap_labels <- c("VD Saprobity 1", "VD Saprobity 2", "VD Saprobity 3", "VD Saprobity 4","VD Saprobity 5", "VD Saprobity Indet", "VD Saprobity Taxa used")
  vdamSaprobity <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(vdamSaprobity) <- vdamsap_labels
  #finds the column
  vdamsap_v <- (taxaInRA[,"vdamsap"])
  #remove the NA
  vdamsap_v[is.na(vdamsap_v)] = 0
  print("Calculating Van Dam saprobity")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamsap_cat1 <- sum(taxaInRA[which(vdamsap_v == 1),sampleNumber])
    vdamsap_cat2 <- sum(taxaInRA[which(vdamsap_v == 2),sampleNumber])
    vdamsap_cat3 <- sum(taxaInRA[which(vdamsap_v == 3),sampleNumber])
    vdamsap_cat4 <- sum(taxaInRA[which(vdamsap_v == 4),sampleNumber])
    vdamsap_cat5 <- sum(taxaInRA[which(vdamsap_v == 5),sampleNumber])
    vdamsap_indet <- 100 - sum(vdamsap_cat1, vdamsap_cat2, vdamsap_cat3, vdamsap_cat4, vdamsap_cat5) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    vdamsaptaxaused <- length(which(vdamsap_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamsaptaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    vdamsap_values <- c(vdamsap_cat1, vdamsap_cat2, vdamsap_cat3, vdamsap_cat4, vdamsap_cat5, vdamsap_indet, vdamsaptaxaused)
    vdamSaprobity[sampleNumber, ] <- vdamsap_values
    if (vandamReports==TRUE & exists("resultsPath")) {write(paste("SAPROBITY- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
    if (vandamReports==TRUE & exists("resultsPath")) {write(Vdamsaptaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
  }
  ##### END VANDAM SAPROBITY

  ##### VANDAM MOISTURE (AERO) STATE
  #creates results dataframe
  vdamaero_labels <- c("VD Aero 1", "VD Aero 2", "VD Aero 3", "VD Aero 4","VD Aero 5", "VD Aero Indet", "VD Aero Taxa used")
  vdamAero <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(vdamAero) <- vdamaero_labels
  #finds the column
  vdamaero_v <- (taxaInRA[,"vdamaero"])
  #remove the NA
  vdamaero_v[is.na(vdamaero_v)] = 0
  print("Calculating Van Dam aero")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamaero_cat1 <- sum(taxaInRA[which(vdamaero_v == 1),sampleNumber])
    vdamaero_cat2 <- sum(taxaInRA[which(vdamaero_v == 2),sampleNumber])
    vdamaero_cat3 <- sum(taxaInRA[which(vdamaero_v == 3),sampleNumber])
    vdamaero_cat4 <- sum(taxaInRA[which(vdamaero_v == 4),sampleNumber])
    vdamaero_cat5 <- sum(taxaInRA[which(vdamaero_v == 5),sampleNumber])
    vdamaero_indet <- 100 - sum(vdamaero_cat1, vdamaero_cat2, vdamaero_cat3, vdamaero_cat4, vdamaero_cat5) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    vdamaerotaxaused <- length(which(vdamaero_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamaerotaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    vdamaero_values <- c(vdamaero_cat1, vdamaero_cat2, vdamaero_cat3, vdamaero_cat4, vdamaero_cat5, vdamaero_indet, vdamaerotaxaused)
    vdamAero[sampleNumber, ] <- vdamaero_values
    if (vandamReports==TRUE & exists("resultsPath")) {write(paste("MOISTURE- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
    if (vandamReports==TRUE & exists("resultsPath")) {write(Vdamaerotaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
  }
  ##### END VANDAM MOISTURE (AERO) STATE

  ##### VANDAM TROPHIC STATE
  #creates results dataframe
  vdamtrop_labels <- c("VD Trophic 1", "VD Trophic 2", "VD Trophic 3", "VD Trophic 4","VD Trophic 5", "VD Trophic Indet", "VD Trophic Taxa used")
  vdamTrophic <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(vdamTrophic) <- vdamtrop_labels
  #finds the column
  vdamtrop_v <- (taxaInRA[,"vdamtrop"])
  #remove the NA
  vdamtrop_v[is.na(vdamtrop_v)] = 0
  print("Calculating Van Dam trophic state")
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    #sum abundances per category
    vdamtrop_cat1 <- sum(taxaInRA[which(vdamtrop_v == 1),sampleNumber])
    vdamtrop_cat2 <- sum(taxaInRA[which(vdamtrop_v == 2),sampleNumber])
    vdamtrop_cat3 <- sum(taxaInRA[which(vdamtrop_v == 3),sampleNumber])
    vdamtrop_cat4 <- sum(taxaInRA[which(vdamtrop_v == 4),sampleNumber])
    vdamtrop_cat5 <- sum(taxaInRA[which(vdamtrop_v == 5),sampleNumber])
    vdamtrop_indet <- 100 - sum(vdamtrop_cat1, vdamtrop_cat2, vdamtrop_cat3, vdamtrop_cat4, vdamtrop_cat5) #calculates indetermined
    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    vdamtroptaxaused <- length(which(vdamtrop_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    Vdamtroptaxaused_taxa <- taxaInRA[which(vdams_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    vdamtrop_values <- c(vdamtrop_cat1, vdamtrop_cat2, vdamtrop_cat3, vdamtrop_cat4, vdamtrop_cat5, vdamtrop_indet, vdamtroptaxaused)
    vdamTrophic[sampleNumber, ] <- vdamtrop_values
    if (vandamReports==TRUE & exists("resultsPath")) {write(paste("TROPHIC STATE- Sample:", colnames(taxaInRA[sampleNumber])),paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
    if (vandamReports==TRUE & exists("resultsPath")) {write(Vdamtroptaxaused_taxa,paste(resultsPath, "\\", vandamtaxafile, sep=""), append=TRUE)}
  }
  ##### END VANDAM TROPHIC STATE

  ### VANDAM RESULTS SUMMARY TABLES
  vandam.results <- data.frame(c(vdamSalinity, vdamNHeterotrophy, vdamOxygen, vdamSaprobity, vdamAero, vdamTrophic))

  if (vandamReports==TRUE & exists("resultsPath")) {print(paste("VanDam detailed reports - File exported as ", paste(resultsPath, "\\", vandamtaxafile, sep="")))}
  ### END VANDAM RESULTS

  return(vandam.results)
}
