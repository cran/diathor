#' Calculate size classes for diatoms
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @description
#' The input for these functions is the resulting dataframe obtained from the diat_loadData() function to calculate size classes for diatoms
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi:10.1016/j.ecolind.2019.105951
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

  sc1 <- sc2 <- sc3 <- sc4 <- sc5 <- Size_class_Indet <- Size_Taxa_used <- NULL

  #gets the column named "species", everything before that is a sample
  lastcol <- which(colnames(taxaInEco)=="species")
  print("Calculating size classes")

  taxaInRA <- taxaInEco
  taxaInRA_samples = taxaInRA[, 1:(lastcol - 1)]
  data.table::setDT(taxaInRA_samples)
  data.table::setnafill(taxaInRA_samples, fill = 0)
  rel_abu  <- apply(taxaInRA_samples, 2, function(x)
    round(x / sum(x) * 100, 2))
  taxaInRA <- cbind(rel_abu, taxaInRA[, lastcol:ncol(taxaInRA)])
  sampleNames <- colnames(taxaInRA[1:(lastcol - 1)])
  data.table::setDT(taxaInRA)
  size_v <- taxaInRA[, which(grepl(pattern = "classe_de_taille", names(taxaInRA))), with = FALSE]
  size_v[is.na(size_v)] = 0

  size.results <- data.table::data.table(
    sc1 <- unlist(taxaInRA[which(size_v == 1), lapply(.SD, sum, na.rm = TRUE), .SDcols = 1:(lastcol -
                                                                                             1)]),
    sc2 <- unlist(taxaInRA[which(size_v == 2), lapply(.SD, sum, na.rm =
                                                       TRUE), .SDcols = 1:(lastcol - 1)]),
    sc3 <- unlist(taxaInRA[which(size_v == 3), lapply(.SD, sum, na.rm =
                                                       TRUE), .SDcols = 1:(lastcol - 1)]),
    sc4 <- unlist(taxaInRA[which(size_v == 4), lapply(.SD, sum, na.rm =
                                                       TRUE), .SDcols = 1:(lastcol - 1)]),
    sc5 <- unlist(taxaInRA[which(size_v == 5), lapply(.SD, sum, na.rm =
                                                       TRUE), .SDcols = 1:(lastcol - 1)])
  )
  size.results[, `Size_class_Indet` := round(100 - (sc1 + sc2 + sc3 + sc4 + sc5), 1)]
  size.results[, `Size_Taxa_used` := numeric(lastcol - 1)]
  size.results[`Size_class_Indet` < 0, `Size_class_Indet` := 0]
  te2 = purrr::map(.x = 1:(lastcol - 1),
                   .f = ~ taxaInRA_samples[, .x, with = F] > 0 & size_v != 0)
  te3 = matrix(unlist(te2), ncol = lastcol - 1, byrow = TRUE) %>% as.data.frame()
  names(te3) = names(taxaInRA_samples)
  size.results[, `Size_Taxa_used` := colSums(te3)]
  names(size.results) = c(
    "Size class 1",
    "Size class 2",
    "Size class 3",
    "Size class 4",
    "Size class 5",
    "Size class Indet",
    "Size Taxa used"
  )
  size.results <- as.data.frame(size.results) #need to convert it to dataframe explicitly to plot
  return(size.results)

}

