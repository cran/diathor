#' Calculate morphological parameters for diatoms
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @param isRelAb Boolean. If set to 'TRUE' it means that your species' data is the relative abundance of each species per site. If FALSE, it means that it the data corresponds to absolute densities. Default = FALSE
#' @description
#' The input for these functions is the resulting dataframe obtained from the diat_loadData() function to calculate morphological parameters
#' The morphological data (size classes, chlorophlasts) is obtained from the 'Diat.Barcode' project. Besides citing DiaThor, the Diat.Barcode project should also be cited if the package is used, as follows:
#' \itemize{
#' \item Rimet F., Gusev E., Kahlert M., Kelly M., Kulikovskiy M., Maltsev Y., Mann D., Pfannkuchen M., Trobajo R., Vasselon V., Zimmermann J., Bouchez A., 2019. Diat.barcode, an open-access curated barcode library for diatoms. Scientific Reports. https://www.nature.com/articles/s41598-019-51500-6
#' }
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
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
#' morphoResults <- diat_morpho(df)
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @export diat_morpho
#'

###### ---------- FUNCTION FOR MORPHOLOGICAL DATA: DONE ---------- ########
#### IN THIS SECTION WE CALCULATE MORPHOLOGICAL VARIABLES
### INPUT: result created in loadData(). If its in RA, biovolume cannot be calculated
### OUTPUTS: #number of chloroplasts #shape of chloroplasts #biovolume

diat_morpho <- function(resultLoad, isRelAb = FALSE){

  # First checks if species data frames exist. If not, loads them from CSV files
  if(missing(resultLoad)) {
    print("Please run the diat_loadData() function first to enter your species data in the correct format")
    #handles cancel button
    if (missing(resultLoad)){
      stop("Calculations cancelled")
    }
  }

  #taxaIn <- resultLoad[[2]] #abundance in absolute
  #taxaInRA <- resultLoad[[1]] #abundance in relative abundance
  sampleNames <- resultLoad[[3]]
  resultsPath <- resultLoad[[4]]
  taxaIn <- resultLoad[[5]]

  #checks thata taxaIn (taxaIn from diat_Load) has at least recognized some species
  if (nrow(taxaIn)==0){
    print("No species were recognized for morphology calculations")
    print("Morphology data will not be available")
    morphoresultTable <- NULL
    return(morphoresultTable)
  }

  #gets the column named "species", everything before that is a sample

  lastcol <- which(colnames(taxaIn)=="species")

  # Convert taxaIn sample data to Relative Abundance data
  if(isRelAb==FALSE){ #pass parameter when in function
    taxaInRA <- taxaIn
    for (i in 1:nrow(taxaIn)){
      for (j in 1:(lastcol-1)){
        if (is.na(taxaIn[i,j])){
          taxaInRA[i,j] <- 0
        } else {
          taxaInRA[i,j] <- (taxaIn[i,j]*100)/sum(taxaIn[,j])
        }
      }
    }
  } else {taxaInRA <- taxaIn}

  ##### Number of Chloroplasts
  #creates an empty dataframe
  numcloroplastos <- c(unique(taxaInRA[,"chloroplast_number"]))
  numcloroplastos.result <- data.frame(matrix(ncol = length(numcloroplastos), nrow = (lastcol-1)))
  colnames(numcloroplastos.result) <- paste("chloroplasts - ", (numcloroplastos), sep ="")

  #PROGRESS BAR
  print("Calculating chloroplast quantity")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix

    numcloroplastos.ab <- NULL
    #sum abundances per category
    for (i in numcloroplastos){
      abc <- sum(taxaInRA[which(taxaInRA$chloroplast_number == i),sampleNumber])
      #round numbers and remove negatives
      if (abc<0){abc <- 0}
      abc <- round(abc, digits=3)
      numcloroplastos.ab <- c(numcloroplastos.ab, abc)
    }

    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    numcloroplastos.result[sampleNumber, ] <- numcloroplastos.ab
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  ##### END Number of Chloroplasts

  ##### Shape of Chloroplasts
  #creates results dataframe
  shpcloroplastos <- c(unique(taxaInRA[,"chloroplast_shape"]))
  shpcloroplastos.result <- data.frame(matrix(ncol = length(shpcloroplastos), nrow = (lastcol-1)))
  colnames(shpcloroplastos.result) <- paste("shape chloroplasts -", (shpcloroplastos), sep ="")

  #PROGRESS BAR
  print("Calculating chloroplast shapes")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)
  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    shpcloroplastos <- c(unique(taxaInRA[,"chloroplast_shape"]))
    shpcloroplastos.ab <- NULL
    #sum abundances per category
    for (i in shpcloroplastos){
      abc <- sum(taxaInRA[which(taxaInRA$chloroplast_shape == i),sampleNumber])
      #round numbers and remove negatives
      if (abc<0){abc <- 0}
      abc <- round(abc, digits=3)
      shpcloroplastos.ab <- c(shpcloroplastos.ab, abc)
    }

    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    shpcloroplastos.result[sampleNumber, ] <- shpcloroplastos.ab
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)
  ##### END Shape of Chloroplasts

  ##### Biovolume

  #creates result dataframe
  biovol.val.result <- data.frame(matrix(ncol = 1, nrow = (lastcol-1)))
  colnames(biovol.val.result) <- "Total Biovolume"
  #PROGRESS BAR
  print("Calculating biovolume")
  pb <- txtProgressBar(min = 1, max = (lastcol-1), style = 3)

  for (sampleNumber in 1:(lastcol-1)){ #for each sample in the matrix
    biovol <- c(unique(taxaInRA[,"biovolume"]))
    biovol[is.na(biovol)] = 0
    taxaInRA[is.na(taxaInRA$biovolume), "biovolume"] = 0
    biovol.val <- NULL
    #sum abundances*biovolume
    biovol.val <- sum((taxaInRA[which(biovol != 0),sampleNumber])*taxaInRA[which(biovol != 0),"biovolume"])
    #labels and exports dataframe with results, in a single row to add the rest of samples as rows
    biovol.val.result[sampleNumber, ] <- biovol.val
    #update progressbar
    setTxtProgressBar(pb, sampleNumber)
  }
  #close progressbar
  close(pb)



  ##### END Biovolume
  morphoresultTable <- NULL
  morphoresultTable <- data.frame(c( if(exists("numcloroplastos.result")){numcloroplastos.result}, if(exists("shpcloroplastos.result")){shpcloroplastos.result}, if(exists("biovol.val.result")){biovol.val.result}))
  rownames(morphoresultTable) <- sampleNames
  morphoresultTable <- list(as.data.frame(numcloroplastos.result), as.data.frame(shpcloroplastos.result), as.data.frame(if(exists("biovol.val.result")){biovol.val.result}))
  names(morphoresultTable) <- c("numcloroplastos.result", "shpcloroplastos.result", "biovol.val.result")
  return(morphoresultTable)
}


