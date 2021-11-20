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
  numcloroplastos.result <- shpcloroplastos.result <- biovol.val.result <- NULL
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
  ### MODIFIED JJ
  if (isRelAb == FALSE) {
    taxaInRA <- taxaIn
    setDT(taxaInRA)
    setnafill(x = taxaInRA[,1:(lastcol-1)], fill = 0)
    rel_abu  = apply(taxaInRA[,1:(lastcol-1)], 2, function(x)
      round(x / sum(x) * 100, 2))
    taxaInRA = cbind(rel_abu, taxaInRA[, lastcol:ncol(taxaInRA)])
  } else {
    taxaInRA <- taxaIn
  }

  #   ##### Number of Chloroplasts / Shape / Biovolume
  taxaInRA_samples = taxaInRA[,1:(lastcol-1)]
  #   #PROGRESS BAR
  print("Calculating chloroplasts and biovolume")
  pb <- txtProgressBar(min = 1, max = 3, style = 3)
  for (i in 1:3){
    lp_var = switch(i,
                    taxaInRA[, "chloroplast_number"],
                    taxaInRA[, "chloroplast_shape"],
                    taxaInRA[, "biovolume"])

    prefix_name = switch(i, "chloroplasts -",
                         "shape chloroplasts -",
                         "Total Biovolume")
    final_name = switch(i,
                        "numcloroplastos.result",
                        "shpcloroplastos.result",
                        "biovol.val.result")
    if (i == 3){
      lp_data = taxaInRA_samples * unlist(lp_var)
      lp_data = lp_data[, lapply(.SD, sum, na.rm = T), .SDcols = 1:(lastcol-1)]
    } else {
      lp_data = taxaInRA_samples[, lapply(.SD, sum, na.rm = T), by = lp_var]
    }
    lp_data = data.table::transpose(lp_data)
    if (i == 3){
      names(lp_data) = prefix_name
    } else {
      names(lp_data) = paste(prefix_name,unlist(lp_data[1,]))
      lp_data = lp_data[-1,]
    }
    #if (i == 1) {
    if (i < 3) { #same procedure for i = 1 and i =2
      lp_data = lp_data[, lapply(.SD, as.numeric)]
      lp_data = lp_data[, lapply(.SD, round, 4)]
    }
    rownames(lp_data) = sampleNames
    assign(x = final_name,
           value = lp_data)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  #Outputs
  names(numcloroplastos.result)
  morphoresultTable <- list(as.data.frame(numcloroplastos.result),
                            as.data.frame(shpcloroplastos.result), as.data.frame(if (exists("biovol.val.result")) {
                              biovol.val.result
                            }))
  names(morphoresultTable) <- c("numcloroplastos.result",
                                "shpcloroplastos.result", "biovol.val.result")
  return(morphoresultTable)




}


