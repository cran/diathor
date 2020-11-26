#' Loads the Data into DiaThor in the correct format
#' @param species_df The data frame with your species data. Species as rows, Samples as columns. If empty, a dialog box will prompt to import a CSV file
#' @param isRelAb Boolean. If set to 'TRUE' it means that your species' data is the relative abundance of each species per site. If FALSE, it means that it the data corresponds to absolute densities. Default = FALSE
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
#' @param resultsPath String. Path for the output data. If empty (default), it will prompt a dialog box to select an output folder
#' @description
#' Loads the CSV or dataframe file, sets the Output folder for the package, and conducts both an exact and an heuristic search of the species' names.
#'
#' The input file for the package is a dataframe or an external CSV file. Species should be listed as rows, with species' names in column 1 (column name should be "species")
#' If the input data contains a column named "acronym", the package will use that column to match species with their ecological values. This is more accurate than the
#' heuristic search of species' names.
#' The other columns (samples) have to contain the abundance of each species (relative or absolute) in each sample.
#' The first row of the file has to contain the headers with the sample names. Remember that a column named "species" is mandatory, containing the species' names
#' If a dataframe is not specified as a parameter (species_df), the package will show a dialog box to search for the CSV file
#' A second dialog box will help set up an Output folder, where all outputs from the package will be exported to (dataframes, CSV files, plots in PDF)
#' The package also downloads and installs a wrapper for the 'Diat.Barcode' project. Besides citing the DiaThor package, the Diat.Barcode project should also be cited, as follows:
#' \itemize{
#' \item Rimet F., Gusev E., Kahlert M., Kelly M., Kulikovskiy M., Maltsev Y., Mann D., Pfannkuchen M., Trobajo R., Vasselon V., Zimmermann J., Bouchez A., 2019. Diat.barcode, an open-access curated barcode library for diatoms. Scientific Reports. https://www.nature.com/articles/s41598-019-51500-6
#' }
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @importFrom stringdist stringdist ain
#' @export diat_loadData

########-------- FUNCTION TO LOAD FILE  --------------#########
#### IN THIS SECTION WE READ THE MATRIX
### INPUT: CSV file with the samples * abundances or relative abundances
### OUTPUTS: a dataframe with the species as matched against the database, with species in RA; Taxaincluded and Taxaexcluded in CSV in the Output
### folder, detailing which taxa were recognized and which were not

diat_loadData <- function(species_df, isRelAb=FALSE, maxDistTaxa=2, resultsPath){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(species_df)) {
    print("Select CSV matrices with your sample data")
    Filters <- matrix(c("Comma Separated Values (CSV)", "*.csv"),
                      1, 2, byrow = TRUE)
    species_df <- as.data.frame(read.csv(file.choose())) #select the data matrix with a dialog box
  }

  #Results folder
  if(missing(resultsPath)) {
    print("Select Results folder")
    resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  }
  if (is.na(resultsPath)){stop("Calculations cancelled, no folder selected")}
  print("Result folder selected")

  #check for duplicate species
  if("species" %in% colnames(species_df)) {
    print("Checking for duplicates")
    testdupl <- species_df[,"species"]
    if (length(testdupl[duplicated(testdupl)])> 0){ #check for duplicates
      print("Duplicate species found in input data")
      testdupl[duplicated(testdupl)] #show duplicates
      stop("Cancelling.Duplicate species found in input data") #abort script
    }
    speciesrow <- species_df[,"species"]
    row.names(species_df) <- speciesrow  #converts the species names to rownames
    species_df<- species_df[ , !(names(species_df) == "species")] #removes the species names column
  } else {
    #check if species names are in the first column
    stop("Calculations cancelled. File should contain a column named 'species' with the species list")
  }
  #Reorder alphabetically
  #species_df <- species_df[order(row.names(species_df)),] #reorder alphabetically


  # Find acronyms or update then
  print("Finding and updating acronyms")
  species_df <- diat_findAcronyms(species_df, maxDistTaxa, resultsPath)

  #removes NA in abundance data
  species_df[is.na(species_df)] <- 0

  ########## LINK WITH DIAT.BARCODE DATABASE
  #internal 'Diat.Barcode' number
  intversion <- "9"
  #get version number of latest 'Diat.Barcode'
  dic <- read.csv("http://www.francoiskeck.fr/work/diatbarcode/dic_version.csv", header = TRUE, stringsAsFactors = FALSE)
  if (exists("dic")){
    #is able to check the version
    version <- dic[dic$Flavor =="original",]
    version <- version$Version[which.max(as.numeric(as.POSIXlt(version$Date, format = "%d-%m-%Y")))]

    #compare both. If updates are needed, attempt them
    if (version == intversion){
      #updates are not needed
      print ("No updates needed for the 'Diat.barcode' database. Proceeding")
      #load("data/dbc_offline.RData") ##takes the internal database
      #dbc <- dbc_offline
      dbc <- diathor::dbc_offline

    } else {
      #updates are needed
      ########--------  Diat.Barcode download attempt. If it fails, tries to use internal database
      ## WARNING: CRAN package does not auto-update the Diat.Barcode database
      print("The diatom database in DiaThor is out of date")
      print("The CRAN version of the package does not auto-update the internal database. But the GitHub version does!")
      print("Using internal database, 'Diat.barcode' v.9.0 published on 14-09-2020")
      dbc <- diathor::dbc_offline

      ###### THIS SECTION IS FOR THE GITHUB PROJECT ONLY
      #print("Attempting to download diat.barcode from website")
      # dbc <- diatbarcode::get_diatbarcode(version = "last") #loads the latest version of diat.barcode
      # if (exists("dbc")){ #it if was able to download the new version, proceed
      #   print("Latest version of Diat.barcode succesfully downloaded. Remember to credit accordingly!")
      # } else { #it if was unable to download the new version, use the internal database
      #   print("Latest version of Diat.barcode cannot be downloaded")
      #   print("Using internal database, Diat.barcode v.8.1 published on 10-06-2020. It might need to be updated")
      #   #load("data/dbc_offline.RData") ##takes the internal database
      #   #dbc <- dbc_offline
      #   dbc <- dbc_offline
      # }
      ###### END OF GITHUB VERSION

    }
  } else {
    print("Latest version of 'Diat.barcode' unknown")
    print("Using internal database, 'Diat.barcode' v.9.0 published on 14-09-2020")
    dbc <- diathor::dbc_offline
  }
  ### Double checks that database got loaded correctly or cancels alltogether
  if (exists("dbc")){
  } else { #if everything fails, cancels
    stop("Database could not be downloaded or retrieved from internal package. Cancelling calculations")
  }
  ########## END LINK WITH DIAT.BARCODE DATABASE

  #Remove duplicate by field "species" in diat.barcode
  dbc2 <- as.data.frame(dbc[!duplicated(dbc[,"species"]),]) #transforms dbc to a dataframe
  ecodata <- dbc2[which(colnames(dbc2)=="species" ):ncol(dbc2)] #keeps only from the "species" column onwards, to keep the ecological data

  ########-------- MATCHES AND BINDS DATA SETS
  # NOTE: if several species match with the same name in the diat.barcode database (multiple strains), it keeps the
  # first occurence. Ecological indices usually match in all strains of the same species
  #fuzzy matching

  #NEW SEARCH
  rownames(species_df) <- trimws(rownames(species_df))
  taxaInSp <- as.data.frame(matrix(nrow=nrow(species_df), ncol = (ncol(species_df)+ncol(ecodata))))
  #species_df section
  taxaInSp[1:nrow(species_df),1:ncol(species_df)] <- species_df[1:nrow(species_df),1:ncol(species_df)] #copies species_df into taxaInSp
  colnames(taxaInSp)[1:ncol(species_df)] <- colnames(species_df) #copies column names
  rownames(taxaInSp) <- rownames(species_df) #copies row names
  taxaInSp$recognizedSp <- NA #creates a new column with the recognized species

  #ecodata_section
  lastcolspecies_df <-  which(colnames(taxaInSp)=="new_species") #gets last column of taxaInSp with species_df data
  colnames(taxaInSp)[(lastcolspecies_df+1):(ncol(taxaInSp)-1)] <- colnames(ecodata) #copies column names
  #PROGRESS BAR
  pb <- txtProgressBar(min = 1, max = nrow(taxaInSp), style = 3)
  for (i in 1:nrow(taxaInSp)){
    searchvectr <- ecodata[stringdist::ain(ecodata[,"species"],row.names(species_df)[i], maxDist=maxDistTaxa, matchNA = FALSE),] #seaches species by species
    if (nrow(searchvectr)==1){ #if it finds only one species, add that
      taxaInSp[i,(lastcolspecies_df + 1):(ncol(taxaInSp)-1)] <- searchvectr
      taxaInSp[i,"recognizedSp"] <- searchvectr$species
    } else if (nrow(searchvectr)>1){ #if it finds multiple species, keeps the one with the lower distance
      searchvectr <- searchvectr[which(stringdist::stringdist(searchvectr$species, row.names(species_df)[i]) == min(stringdist::stringdist(searchvectr$species, row.names(species_df)[i]))),]
      if (nrow(searchvectr) > 1) { #still finds more than one with the same lower distance, keeps the first
        searchvectr <- searchvectr[1,]
      }
      taxaInSp[i,(lastcolspecies_df + 1):(ncol(taxaInSp)-1)] <- searchvectr
      taxaInSp[i,"recognizedSp"] <- searchvectr$species

    } else if (nrow(searchvectr)==0){ #species not found at all
      taxaInSp[i,"recognizedSp"] <- "Not found"
    }
    #update progressbar
    setTxtProgressBar(pb, i)
  }
  #close progressbar
  close(pb)
  #END NEW SEARCH
  taxaInEco <- taxaInSp
  taxaIn <- species_df #dataframe to be exported for indices
  taxaIncluded <- as.data.frame(rownames(taxaInEco)[which(taxaInEco$recognizedSp != "Not found")])
  taxaExcluded <- as.data.frame(rownames(taxaInEco)[which(taxaInEco$recognizedSp == "Not found")])

  #remove the acronym col and move it to the back
  acronym_col <- taxaIn$acronym
  acronym_col2 <- taxaInEco$acronym
  acronym_col <- replace(acronym_col, acronym_col=="0", NA)
  acronym_col2 <- replace(acronym_col2, acronym_col2=="0", NA)
  taxaIn<- taxaIn[ , !(names(taxaIn) == "acronym")] #removes the acronym names column
  taxaInEco<- taxaInEco[ , !(names(taxaInEco) == "acronym")] #removes the acronym names column
  taxaInEco$acronym <- acronym_col2
  taxaIn$acronym <- acronym_col

  #remove the updated species col and move it to the back in both taxaIn
  newspecies_col <- taxaIn$new_species
  newspecies_col <- replace(newspecies_col, newspecies_col=="0", NA)
  taxaIn<- taxaIn[ , !(names(taxaIn) == "new_species")] #removes the newspecies_col
  taxaIn$new_species <- newspecies_col

  newspecies_col2 <- taxaInEco$new_species
  newspecies_col2 <- replace(newspecies_col2, newspecies_col2=="0", NA)
  taxaInEco<- taxaInEco[ , !(names(taxaInEco) == "new_species")] #removes the newspecies_col
  taxaInEco$new_species <- newspecies_col2

  sampleNames <- colnames(taxaIn) #Saves the names of the samples
  removeelem <- c("acronym","new_species","species") #Removes columns not samples
  sampleNames <- sampleNames[!(sampleNames %in% removeelem)]


  if (nrow(taxaIncluded) == 0) {
    print("No taxa were recognized for morphology analysis. Check taxonomy instructions for the package")
  }
  #Exports included taxa in morphology analyses
  colnames(taxaIncluded) <- "Eco/Morpho"
  write.csv(taxaIncluded, paste(resultsPath,"\\Taxa included.csv", sep=""))
  print(paste("Number of taxa recognized for morphology:", nrow(taxaIncluded), "-- Detailed list in 'Taxa included.csv'"))

  #also makes a matrix for all the taxa left out, for the user to review
  colnames(taxaExcluded) <- "Eco/Morpho"
  write.csv(taxaExcluded, paste(resultsPath,"\\Taxa excluded.csv", sep=""))
  print(paste("Number of taxa excluded for morphology:", nrow(taxaExcluded), "-- Detailed list in 'Taxa excluded.csv'" ))

  #Has to clean the TaxaInEco for those species that were not found
  taxaInEco <- taxaInEco[which(taxaInEco$recognizedSp != "Not found"),]

  #creates a blank precision matrix for the indices
  precisionmatrix <- as.data.frame(sampleNames)
  names(precisionmatrix)[names(precisionmatrix)=="sampleNames"] <- "Sample"

  write.csv(precisionmatrix, paste(resultsPath,"\\Precision.csv", sep=""))

  #gets the column named "acronym" everything before that column should be a sample with abundance data
  lastcol = which(colnames(taxaIn)=="acronym")
  #
  # #gets sample names
  # sampleNames <- colnames(taxaIn[1:(lastcol-1)])

  #CREATES A RELATIVE ABUNDANCE MATRIX AS WELL FOR THOSE INDICES THAT USE IT
  #Convert taxaIn sample data to Relative Abundance data
  if(isRelAb==FALSE){ #pass parameter when in function
    taxaInRA <- taxaIn
    print("Converting species' densities to relative abundance")
    for (i in 1:nrow(taxaIn)){
      for (j in 1:(lastcol-1)){
        if (is.na(taxaIn[i,j])){
          taxaInRA[i,j] <- 0
        } else {
          taxaInRA[i,j] <- (taxaIn[i,j]*100)/sum(taxaIn[,j])
        }
      }
    }
  } else {
    taxaInRA <- taxaIn
  }

  #CREATES THE EXPORT PRODUCTS
  resultList <- list(as.data.frame(taxaInRA), as.data.frame(taxaIn), sampleNames, resultsPath, taxaInEco)
  names(resultList) <- c("taxaInRA", "taxaIn", "sampleNames", "resultsPath", "taxaInEco")
  print("Data loaded")
  return(resultList)

}
