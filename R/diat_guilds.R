#' Calculate ecological guilds for diatoms
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for these functions is the resulting dataframe obtained from the diat_loadData() function, to calculate the ecological guilds for the diatoms
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
#' }
#' Guild classification is obtained from:
#' \itemize{
#' \item Rimet F. & Bouchez A., 2012. Life-forms, cell-sizes and ecological guilds of diatoms in European rivers. Knowledge and management of aquatic ecosystems, 406: 1-14. DOI:10.1051/kmae/2012018
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
#' guildsResults <- diat_guilds(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_guilds
#'

###### ---------- FUNCTION FOR ECOLOGICAL GUILDS: DONE   ---------- ########
#### IN THIS SECTION WE CALCULATE ECOLOGICAL GUILDS ACCORDING TO PASSY (2017)
### INPUT: taxaInRA created in loadData()
### OUTPUTS: dataframe with ecological guilds  per sample

diat_guilds <- function(resultLoad){
  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaInEco <- resultLoad[[5]]


  #checks thata taxaInEco (taxaInEco from diat_Load) has at least recognized some species
  if (nrow(taxaInEco)==0){
    print("No species were recognized for guild calculations")
    print("Guild data will not be available")
    guilds.results <- NULL
    return(guilds.results)
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

  #removes NA from taxaIn
  taxaInRA[is.na(taxaInRA)] <- 0

  #creates results dataframe
  guild_labels <- c("Guild: HP", "Guild: LP", "Guild: Mot", "Guild: Plank", "Guild: Indet", "Guilds Taxa used")
  guilds.results <- data.frame(matrix(ncol = 6, nrow = (lastcol-1)))
  colnames(guilds.results) <- guild_labels

  #PROGRESS BAR
  print("Calculating ecological guilds")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix

    guild_HP <- taxaInRA[,startsWith(colnames(taxaInRA), "high_profile_guild")]
    guild_LP <- taxaInRA[,startsWith(colnames(taxaInRA), "low_profile_guild")]
    guild_Mot <- taxaInRA[,startsWith(colnames(taxaInRA), "motile_guild")]
    guild_Plank <- taxaInRA[,startsWith(colnames(taxaInRA), "euplanctonic_guild")]
    #remove the NA
    guild_HP[is.na(guild_HP)] = 0
    guild_LP[is.na(guild_LP)] = 0
    guild_Mot[is.na(guild_Mot)] = 0
    guild_Plank[is.na(guild_Plank)] = 0

    #sum abundances per guild
    guild_HP_ab <- sum(taxaInRA[which(guild_HP == 1),sampleNumber])
    guild_LP_ab <- sum(taxaInRA[which(guild_LP == 1),sampleNumber])
    guild_Mot_ab <- sum(taxaInRA[which(guild_Mot == 1),sampleNumber])
    guild_Plank_ab <- sum(taxaInRA[which(guild_Plank == 1),sampleNumber])
    guild_indet <- 100 - sum(guild_HP_ab, guild_LP_ab, guild_Mot_ab, guild_Plank_ab) #calculates indetermined


    #round numbers and remove negatives
    if (guild_indet<0){guild_indet <- 0}
    guild_HP_ab <- round(guild_HP_ab, digits=3)
    guild_LP_ab <- round(guild_LP_ab, digits=3)
    guild_Mot_ab <- round(guild_Mot_ab, digits=3)
    guild_Plank_ab <- round(guild_Plank_ab, digits=3)
    guild_indet <- round(guild_indet, digits=3)

    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    guildtaxaused <- length(which(guild_HP == 1 & taxaInRA[,sampleNumber] > 0))
    guildtaxaused <- guildtaxaused + length(which(guild_LP == 1 & taxaInRA[,sampleNumber] > 0))
    guildtaxaused <- guildtaxaused + length(which(guild_Mot == 1 & taxaInRA[,sampleNumber] > 0))
    guildtaxaused <- guildtaxaused + length(which(guild_Plank == 1 & taxaInRA[,sampleNumber] > 0))

    #which taxa were used? to export
    guildtaxaused_taxa <- taxaInRA[which(guild_HP == 1 & taxaInRA[,sampleNumber] > 0),"species"]
    guildtaxaused_taxa <- c(guildtaxaused_taxa, taxaInRA[which(guild_LP == 1 & taxaInRA[,sampleNumber] > 0),"species"])
    guildtaxaused_taxa <- c(guildtaxaused_taxa, taxaInRA[which(guild_Mot == 1 & taxaInRA[,sampleNumber] > 0),"species"])
    guildtaxaused_taxa <- c(guildtaxaused_taxa, taxaInRA[which(guild_Plank == 1 & taxaInRA[,sampleNumber] > 0),"species"])

    #labels and exports dataframe with results, in a single row to add the rest of samples as rows

    guild_values <- c(guild_HP_ab, guild_LP_ab, guild_Mot_ab, guild_Plank_ab, guild_indet, guildtaxaused)
    guilds.results[sampleNumber, ] <- guild_values
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)

  ##### END GUILDS
  return(guilds.results)

}

