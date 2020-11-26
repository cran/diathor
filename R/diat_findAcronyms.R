#' Search species' names or acronyms for diatoms
#' @param species_df The data frame with your species data. Species as rows, Samples as columns. If empty, a dialog box will prompt to import a CSV file
#' @param resultsPath String. Path to the output folder. If none specified (default), a dialog box will prompt to select it
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
#' @description
#' This function conducts both an exact and an heuristic search of the species' names and tries to convert it to its acronym in the internal database. If acronyms are already present in the input data, it attempts to update them to the latest taxonomy.
#' The input file for the package is a dataframe or an external CSV file. Species should be listed as rows, with species' names in column 1 (column name should be "species")
#' If the input data contains a column named "acronym", the package will use that column to match species with their ecological values. This is more accurate than the
#' heuristic search of species' names.
#' The other columns (samples) have to contain the abundance of each species (relative or absolute) in each sample.
#' The first row of the file has to contain the headers with the sample names. Remember that a column named "species" is mandatory, containing the species' names
#' If a dataframe is not specified as a parameter (species_df), the package will show a dialog box to search for the CSV file
#' A second dialog box will help set up an Output folder, where all outputs from the package will be exported to (dataframes, CSV files, plots in PDF)
#'
#'
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @importFrom stringdist amatch
#'
#' @export diat_findAcronyms

########-------- FUNCTION TO FIND ACRONYMS FROM DATAFRAME  --------------#########
### INPUT: CSV file with the list of species in the "species" column, or dataframe with the same column
### OUTPUT: dataframe with an additional "acronym" column. NA = acronym not found

diat_findAcronyms <- function(species_df, maxDistTaxa=2, resultsPath){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(species_df)) {
    print("Select CSV matrices with your sample data")
    Filters <- matrix(c("Comma Separated Values (CSV)", "*.csv"),
                      1, 2, byrow = TRUE)
    species_df <- as.data.frame(read.csv(file.choose())) #select the data matrix with a dialog box
    #handles cancel button
    if (missing(species_df)){
      stop("Calculations cancelled")
    }
  }

  #Results fomder
  if(missing(resultsPath)) {
    print("Select Results folder")
    resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  }
  if (is.na(resultsPath)){stop("Calculations cancelled, no folder selected")}

  acronymDB <- as.data.frame(diathor::acronyms) #loads taxonomy
  #acronymDB <- read.csv("Indices/acronyms.csv") #uses the external csv file


  #trims acronyms and converts them to characters, so it doesnt recognize them as factors (darn R)
  acronymDB$species <- as.character(acronymDB$species)
  acronymDB$species <- trimws(acronymDB$species)
  acronymDB$acronym <- as.character(acronymDB$acronym)
  acronymDB$acronym <- trimws(acronymDB$acronym)
  acronymDB$new_acronym <- as.character(acronymDB$new_acronym)
  acronymDB$new_acronym <- trimws(acronymDB$new_acronym)


  #FINDS ACRONYMS FROM SPECIES OR PROCEEDS TO ACRONYM UPDATE
  if("acronym" %in% colnames(species_df)) {
    #there is an acronym, updates
    print("Acronyms column found in file. Updating acronyms")
    species_df$acronym <- as.character(species_df$acronym)

  }else{
    #there is no acronym column, creates one
    species_df$acronym <- NA

    print ("Searching for acronyms in exact match, please wait...")
    #exact match by species name (fast)
    species_df$acronym <- acronymDB$acronym[match(trimws(rownames(species_df)), acronymDB$species)]
    #species_df_pre <- species_df
    #species_df <- species_df_pre
    #with the unmatched exact, try heuristic
    print ("Searching for acronyms in heuristic match, please wait...")
    #PROGRESS BAR FOR HEURISTIC SEARCH OF ACRONYMS
    ntot <- nrow(species_df)
    pb <- txtProgressBar(min = 1, max = ntot, style = 3)
    for(i in 1:ntot) {
      if (is.na(species_df$acronym[i])){
        species_df$acronym[i] <- acronymDB$acronym[stringdist::amatch(trimws(rownames(species_df[i,])), trimws(acronymDB$species), maxDist=maxDistTaxa, matchNA = FALSE)]
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }

  #print(paste("Updated acronyms: ", updatedacronyms))
  print(paste("Acronyms were found for ", nrow(species_df) - sum(is.na(species_df$acronym)) , " species out of ", length(species_df$acronym), sep=""))

  #original species and acronym list
  speciesList <- trimws(rownames(species_df))
  acronymList <- trimws(species_df$acronym)
  #builds column with new name for species
  species_df$new_species <- NA

  #UPDATE ACRONYMS
  updatedacronyms <- 0
  species_df$acronym <- as.character(species_df$acronym)

  #keepupdating <- T
  #while (keepupdating == TRUE){
  updatedacronyms <- 0
  for(i in 1:nrow(species_df)) {
    if (!is.na(species_df$acronym[i])){ #there is an acronym for this species, try updating
      updatedacronym <- acronymDB$new_acronym[match(species_df$acronym[i], acronymDB$acronym)] #get updated acronym
      updatedspecies <- acronymDB$new_species[match(species_df$acronym[i], acronymDB$acronym)] #get updated species
      if (is.na(updatedacronym) | updatedacronym == ""){
        #print(paste("i=", i, " already exists"))

      } else{
        #species_df$new_acronym[i] <- updatedacronym #update the acronym
        species_df$new_species[i] <- updatedspecies #update the species
        updatedacronyms <- updatedacronyms + 1
      }
    }
  }
  if (updatedacronyms > 0){
    print(paste(updatedacronyms, "updated acronyms were found and exported"))

  }

  #if (updatedacronyms==0){keepupdating <- F}
  #}



  #taxa recognized by acronym
  print("Exporting list of species with the acronyms found by heuristic search. NA= not found")
  outacronyms <- cbind(speciesList, species_df$new_species, acronymList, species_df$acronym)
  # outacronyms <- cbind(speciesList, acronymList, species_df$acronym)

  colnames(outacronyms) <- c("Original species","Updated species" ,"Original Acronym", "Updated acronym")
  write.csv(outacronyms, paste(resultsPath,"\\Recognized_acronyms.csv", sep=""))
  return(species_df)

}

