#' Loads the Data into DiaThor in the correct format
#' @param species_df The data frame with your species data. Species as rows, Samples as columns. If empty, a dialog box will prompt to import a CSV file
#' @param isRelAb Boolean. If set to 'TRUE' it means that your species' data is the relative abundance of each species per site. If FALSE, it means that it the data corresponds to absolute densities. Default = FALSE
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
#' @param resultsPath String. Path for the output data. If empty (default), it will prompt a dialog box to select an output folder
#' @description
#' Loads the CSV or dataframe file, sets the Output folder for the package, and conducts both an exact and an heuristic search of the species' names.
#'
#' The input file for the package is a dataframe or an external CSV file. Species should be listed as rows, with species' names in column 1 (column name should be "species")
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
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi:10.1016/j.ecolind.2019.105951
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @importFrom stringdist stringdist ain
#' @export diat_loadData

########-------- FUNCTION TO LOAD FILE  --------------#########
#### This function inputs the data and runs the species recognition
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

  #Choose a Results folder
  if(missing(resultsPath)) {
    print("Select Results folder")
    resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  }
  if (is.na(resultsPath)){stop("Calculations cancelled, no folder selected")}
  print("Result folder selected")

  #Check for duplicate species
  if("species" %in% colnames(species_df)) {
    print("Checking for duplicates")
    testdupl <- species_df[,"species"]
    if (length(testdupl[duplicated(testdupl)])> 0){ #check for duplicates
      #duplicates found
      print("Duplicate species found in input data")
      testdupl[duplicated(testdupl)] #show duplicates
      stop("Cancelling.Duplicate species found in input data") #abort script
    } else {
      print("No duplicate species found in input data")
    }
    row.names(species_df) <- species_df[,"species"]  #converts the species names to rownames
    species_df<- species_df[ , !(names(species_df) == "species")] #removes the species names column
  } else {
    #check if species names are in the first column
    stop("Calculations cancelled. File should contain a column named 'species' with the species list")
  }

  #Find and report any empty samples
  colvector <- c()
  for (i in 1:ncol(species_df)){
    if (is.numeric(species_df[,i])){
      if (colSums(species_df[i], na.rm=TRUE) == 0){
        colvector <- c(colvector, i)
      }
    }
  }
  if (length(colvector) != 0){
    if (!is.na(colvector)){
      print(paste("Empty sample removed, column #", colvector))
    }
    species_df <- species_df[,-colvector] #entry dataframe without empty samples

  }

  #If the acronym column exists, removes it (for compatibility issues)
  acrocol <- which(colnames(species_df)=="acronym") #finds the acronym column
  if (length(acrocol) != 0) {
    species_df1 <- species_df[,-acrocol] #removes the acronym column if it exists
    sampleNames <- colnames(species_df1) #Saves the names of the samples

  } else {
    sampleNames <- colnames(species_df) #Saves the names of the samples
  }

  #SINCE THE ACRONYM COLUMN WAS REMOVED, WE HAVE TO CREATE A new_species COLUMN

  species_df$new_species <- NA
  species_df$new_species <- rownames(species_df)

  #removes NA in abundance data
  species_df[is.na(species_df)] <- 0


  ########## LINK WITH DIAT.BARCODE DATABASE (v.0.0.8)
  getDiatBarcode <- diathor::diat_getDiatBarcode() #function that gets the Diat.Barcode database
  ecodata <- as.data.frame(getDiatBarcode[1]) #ecodata
  taxaList <- as.data.frame(getDiatBarcode[2]) #taxaList: diat_getDiatBarcode uses the taxaList() function to build a single list with all taxa. This is the result

  #Trim blank spaces from input data
  rownames(species_df) <- trimws(rownames(species_df))

  ###### START NAME CHECK

  #PROGRESS BAR - SEARCH SPECIES FOR POSSIBLE MISSPELLS
  pb <- txtProgressBar(min = 1, max = nrow(species_df), style = 3)
  #macthing function, uses stringdist package
  for (i in 1:nrow(species_df)){
    #get the species name
    spname <- row.names(species_df)[i]
    #remove double spaces and trim blank characters
    spname <- gsub("(?<=[\\s])\\s*|^\\s+|\\s+$", "", spname, perl = TRUE)  # Remove spaces

    #Species that have "var." or "fo.", the period should be followed by a space. If not, add it
    getPeriods <- lapply(strsplit(spname, ''), function(x) which(x == '.'))
    #if there is a period
    if (length(getPeriods[[1]]) > 0){
      for (j in 1:length(getPeriods)){
        #get index
        periodIndex <- getPeriods[[j]]
        if (substr(spname, periodIndex +1, periodIndex+1)!= " "){
          #if no space after the period, add it
          spname <- paste(substring(spname, c(1,periodIndex+1), c(periodIndex-1,nchar(spname))), collapse=". ")
        }
      }
    }

    #then correct the species name for all future uses
    row.names(species_df)[i] <- spname

    #searches heuristically for the species name agains the full list of taxaList (all taxa in all internal databases)
    searchvectr <- taxaList[stringdist::ain(tolower(taxaList[,"species"]),tolower(spname), maxDist=maxDistTaxa, matchNA = FALSE),] #seaches species by species
    if (length(searchvectr)==0){ #if species not found in any db
      print(paste("Taxon not found in any internal database:", spname))
    }
    #update progressbar
    setTxtProgressBar(pb, i)
  }
  #close progressbar
  close(pb)

  ###### END NAME CHECK

  ########-------- MATCHES AND BINDS DATA SETS
  # NOTE: if several species match with the same name in the diat.barcode database (multiple strains), it keeps the
  # first occurence. Ecological indices usually match in all strains of the same species
  #fuzzy matching

  #SEARCH
  #creates empty data frame taxaInSp
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
  #macthing function, uses stringdist package
  for (i in 1:nrow(taxaInSp)){
    searchvectr <- ecodata[stringdist::ain(ecodata[,"species"],row.names(species_df)[i], maxDist=maxDistTaxa, matchNA = FALSE),] #seaches species by species
    if (nrow(searchvectr)==1){ #if it finds only one species, add that
      taxaInSp[i,(lastcolspecies_df + 1):(ncol(taxaInSp)-1)] <- searchvectr
      taxaInSp[i,"recognizedSp"] <- searchvectr$species
    } else if (nrow(searchvectr)>1){ #if it finds multiple species, keeps the one with the lower distance
      searchvectr <- searchvectr[which(stringdist::stringdist(searchvectr$species, row.names(species_df)[i]) == min(stringdist::stringdist(searchvectr$species, row.names(species_df)[i]))),]
      if (nrow(searchvectr) > 1) { #still finds more than one with the same lower distance, creates a majority consensus for each column
        consensus_row <- matrix(nrow = 1, ncol = ncol(searchvectr))
        for (j in 1:ncol(searchvectr)){
          consensus_row[1,j] <- names(which.max(table(searchvectr[,j])))
        }
        searchvectr_names <- colnames(searchvectr)
        searchvectr <- as.data.frame(consensus_row)
        colnames(searchvectr) <- searchvectr_names
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

  #remove the updated species col  and move it to the back in both taxaIn
  newspecies_col <- taxaIn$new_species
  newspecies_col <- replace(newspecies_col, newspecies_col=="0", NA)
  taxaIn<- taxaIn[ , !(names(taxaIn) == "new_species")] #removes the newspecies_col
  taxaIn$new_species <- newspecies_col
  newspecies_col2 <- taxaInEco$new_species
  newspecies_col2 <- replace(newspecies_col2, newspecies_col2=="0", NA)
  taxaInEco<- taxaInEco[ , !(names(taxaInEco) == "new_species")] #removes the newspecies_col
  taxaInEco$new_species <- newspecies_col2

  #If the acronym column exists in taxaIn, moves it to the back
  acrocol <- which(colnames(taxaIn)=="acronym") #finds the acronym column
  if (length(acrocol) != 0) {
    acrocolv <- taxaIn[,acrocol] #saves the vector
    taxaIn <- taxaIn[,-acrocol] #removes the acronym column if it exists
    taxaIn$acronym <- acrocolv
  }

  #If the acronym column exists in taxaInEco, moves it to the back
  acrocol <- which(colnames(taxaInEco)=="acronym") #finds the acronym column
  if (length(acrocol) != 0) {
    acrocolv <- taxaInEco[,acrocol] #saves the vector
    taxaInEco <- taxaInEco[,-acrocol] #removes the acronym column if it exists
    taxaInEco$acronym <- acrocolv
  }

  removeelem <- c("species") #Removes columns not samples
  sampleNames <- sampleNames[!(sampleNames %in% removeelem)]


  if (nrow(taxaIncluded) == 0) {
    print("No taxa were recognized for morphology analysis. Check taxonomy instructions for the package")
  }
  #Exports included taxa in morphology analyses
  colnames(taxaIncluded) <- "Eco/Morpho"
  write.csv(taxaIncluded, file.path(resultsPath,"Taxa included.csv"))
  print(paste("Number of taxa recognized for morphology:", nrow(taxaIncluded), "-- Detailed list in 'Taxa included.csv'"))

  #also makes a matrix for all the taxa left out, for the user to review
  colnames(taxaExcluded) <- "Eco/Morpho"
  write.csv(taxaExcluded, file.path(resultsPath,"Taxa excluded.csv"))
  print(paste("Number of taxa excluded for morphology:", nrow(taxaExcluded), "-- Detailed list in 'Taxa excluded.csv'" ))

  #Has to clean the TaxaInEco for those species that were not found
  taxaInEco <- taxaInEco[which(taxaInEco$recognizedSp != "Not found"),]

  #creates a blank precision matrix for the indices
  precisionmatrix <- as.data.frame(sampleNames)
  names(precisionmatrix)[names(precisionmatrix)=="sampleNames"] <- "Sample"

  write.csv(precisionmatrix, file.path(resultsPath, "num_taxa.csv"))

  #gets the column named "new_species" everything before that column should be a sample with abundance data
  lastcol = which(colnames(taxaIn)=="new_species")


  #CREATES A RELATIVE ABUNDANCE MATRIX AS WELL FOR THOSE INDICES THAT USE IT
  #Convert taxaIn sample data to Relative Abundance data
  if (isRelAb == FALSE) {
    taxaInRA <- taxaIn

    print("Converting species' densities to relative abundance")
    rel_abu  = apply(taxaInRA[,1:(lastcol-1)], 2, function(x)
      round(x / sum(x) * 100, 2))
    taxaInRA <- as.data.frame(cbind(rel_abu, taxaInRA[, lastcol:ncol(taxaInRA)]))
    names(taxaInRA)[lastcol] <- "new_species"
  } else {
    taxaInRA <- taxaIn
  }
  #converts to double taxaInRA
  lastcol <- which(colnames(taxaInRA)=="new_species")
  for (i in 1:(lastcol-1)){
    taxaInRA[,i] <- as.double(taxaInRA[,i])
  }


  #CREATES THE EXPORT PRODUCTS
  resultList <- list(as.data.frame(taxaInRA), as.data.frame(taxaIn), sampleNames, resultsPath, taxaInEco)
  names(resultList) <- c("taxaInRA", "taxaIn", "sampleNames", "resultsPath", "taxaInEco")
  print("Data loaded")
  return(resultList)

}
