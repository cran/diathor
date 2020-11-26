#' Calculate size classes for diatoms
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for these functions is the resulting dataframe obtained from the diat_loadData() function to calculate size classes for diatoms
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
#' }
#' Size class classification is obtained from:
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
#' sizeResults <- diat_size(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_size
#'

###### ---------- FUNCTION FOR SIZE CLASSES: DONE  ---------- ########
#### IN THIS SECTION WE CALCULATE SIZE CLASSES ACCORDING TO....
### INPUT: resultLoad created in loadData()
### OUTPUTS: dataframe with % of size classes

diat_size <- function(resultLoad){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }
  taxaInEco <-  resultLoad[[5]]

  #checks thata taxaInEco (taxaInEco from diat_Load) has at least recognized some species
  if (nrow(taxaInEco)==0){
    print("No species were recognized for size class calculations")
    print("Size class data will not be available")
    size.results <- NULL
    return(size.results)
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



  #gets sample names
  sampleNames <- colnames(taxaInRA[1:(lastcol-1)])

  ##### SIZE CLASSES
  #creates results dataframe
  size_labels <- c("Size class 1", "Size class 2", "Size class 3", "Size class 4", "Size class 5", "Size class Indet", "Size Taxa used")
  size.results <- data.frame(matrix(ncol = 7, nrow = (lastcol-1)))
  colnames(size.results) <- size_labels

  #finds the size class column
  size_v <-  taxaInRA[,startsWith(colnames(taxaInRA), "classe_de_taille")]
  #remove the NA
  size_v[is.na(size_v)] = 0

  #PROGRESS BAR
  print("Calculating size classes")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix

    #sum abundances per category
    size_cat1 <- sum(taxaInRA[which(size_v == 1),sampleNumber])
    size_cat2 <- sum(taxaInRA[which(size_v == 2),sampleNumber])
    size_cat3 <- sum(taxaInRA[which(size_v == 3),sampleNumber])
    size_cat4 <- sum(taxaInRA[which(size_v == 4),sampleNumber])
    size_cat5 <- sum(taxaInRA[which(size_v == 5),sampleNumber])
    size_indet <- 100 - sum(size_cat1, size_cat2, size_cat3, size_cat4, size_cat5) #calculates indetermined
    #round numbers and remove negatives
    if (size_indet<0){size_indet <- 0}
    size_cat1 <- round(size_cat1, digits=3)
    size_cat2 <- round(size_cat2, digits=3)
    size_cat3 <- round(size_cat3, digits=3)
    size_cat4 <- round(size_cat4, digits=3)
    size_cat5 <- round(size_cat5, digits=3)
    size_indet <- round(size_indet, digits=3)

    #how many taxa will be used to calculate? Taxa that have a valid indicator value and abundance > 0
    sizetaxaused <- length(which(size_v != 0 & taxaInRA[,sampleNumber] > 0))
    #which taxa were used? to export
    sizetaxaused_taxa <- taxaInRA[which(size_v != 0 & taxaInRA[,sampleNumber] > 0),"species"]
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    size_values <- c(size_cat1, size_cat2, size_cat3, size_cat4, size_cat5, size_indet, sizetaxaused)
    size.results[sampleNumber, ] <- size_values
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)

  }
  #close progressbar
  close(pb)
  ##### END SIZE CLASSES
  return(size.results)
}

