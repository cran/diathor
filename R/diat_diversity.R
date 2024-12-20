#' Calculate diversity parameters for diatoms using the vegan package
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for these functions is the resulting dataframe obtained from the diat_loadData() function, to calculate diversity data using the vegan package
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi:10.1016/j.ecolind.2019.105951
#' }
#' Diversity index (Shannons H') is calculated using the vegan package, following:
#' \itemize{
#' \item Shannon, C. E., and Weaver, W. (1949). ‘The Mathematical Theory of Communication.’ (University of Illinios Press: Urbana, IL, USA.)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @importFrom vegan specnumber diversity
#' @export diat_diversity
#'

###### ---------- FUNCTION FOR DIVERSITY INDICES: DONE   ---------- ########
#### IN THIS SECTION WE CALCULATE BASIC ECOLOGICAL INDICES WITH VEGAN
### INPUT: resultLoad Cannot be in Relative Abundance
### OUTPUTS: dataframe with diversity indices per sample

diat_diversity <- function(resultLoad){


  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  taxaIn <- resultLoad[[2]]

  #removes NA from taxaIn
  taxaIn[is.na(taxaIn)] <- 0

  #gets the column named "new_species", everything before that is a sample
  lastcol <- which(colnames(taxaIn)=="new_species")

  #remove NA to 0
  taxaIn[is.na(taxaIn)] <- 0

  sampleNames <- colnames(taxaIn[1:(lastcol-1)])
  diversityIndices <- data.frame(matrix(ncol = (lastcol-1), nrow = 3))

  #PROGRESS BAR
  print("Calculating diversity indices")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  #samples are columns

  for (i in 1:(lastcol-1)){

    richness <-vegan::specnumber(taxaIn[,i])
    shannon <- vegan::diversity(taxaIn[,i])
    evenness <- shannon/log(vegan::specnumber(taxaIn[,i]))
    diversityIndices[1,i] <- richness
    diversityIndices[2,i] <- shannon
    diversityIndices[3,i] <- evenness
    #update progressbar
    setTxtProgressBar(pb, i)
  }
  #close progressbar
  close(pb)
  #RESULTS
  diversity.results<-as.data.frame(t(diversityIndices)) #transposes the diversity matrix to plot
  rownames(diversity.results) <- colnames(taxaIn[1:(lastcol-1)])
  colnames(diversity.results) <- c("Richness", "Shannon H", "Evenness")

  return(diversity.results)
}


