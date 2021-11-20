#' Calculates ecological information for diatoms based on the Van Dam classification
#' @param resultLoad The resulting list obtained from the diat_loadData() function
#' @param vandamReports Boolean. If set to 'TRUE' the detailed reports for the Van Dam classifications will be reported in the Output. Default = TRUE
#' @description
#' The input for these functions is the resulting dataframe obtained from the diat_loadData() function, to calculate ecological information for diatoms based on the Van Dam classification
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi:10.1016/j.ecolind.2019.105951
#' }
#' Van Dam classification is obtained form:
#' \itemize{
#' \item Van Dam, H., Mertens, A., & Sinkeldam, J. (1994). A coded checklist and ecological indicator values of freshwater diatoms from the Netherlands. Netherland Journal of Aquatic Ecology, 28(1), 117-133.
#' }
#' @examples
#' \dontrun{
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
  v1 <- v2 <- vd1 <- vd2 <- vd3 <- vd4 <- vd5 <- vd6 <- vd7 <-NULL
  #gets the column named "species", everything before that is a sample
  lastcol <- which(colnames(taxaInEco)=="species")
  taxaInRA <- taxaInEco
  setDT(taxaInRA)
  taxaInRa_samples = taxaInRA[,1:(lastcol-1)]
  setnafill(taxaInRa_samples, fill=0)
  rel_abu  = apply(taxaInRa_samples, 2, function(x)
    round(x / sum(x) * 100, 2))
  taxaInRA = cbind(rel_abu, taxaInRA[, lastcol:ncol(taxaInRA)])
  taxaInRA[is.na(taxaInRA)] <- 0
  resultsPath <- resultLoad[[4]]
  if (vandamReports == TRUE & is.na(resultsPath) == TRUE) {
    print("Select Results folder for the detailed VanDam reports")
    resultsPath <- choose.dir(default = "", caption = "Select folder for your Results")
  }
  if (vandamReports == TRUE & is.na(resultsPath) == FALSE) {
    print("Exporting detailed reports for VanDam ecological preferences")
    print(resultsPath)
    vandamtaxafile = paste("VanDam Taxa used.csv", sep = "")
    write("TAXA USED FOR EACH ECOLOGICAL VARIABLE USING VANDAM's CLASSIFICATION",
          file = file.path(resultsPath,vandamtaxafile))
    write("These taxa were included because: ", file = file.path(resultsPath,vandamtaxafile), append = TRUE)
    write("a) they had a reliable classification in VanDam's classification system", file = file.path(resultsPath,vandamtaxafile), append = TRUE)
    write("b) they had a relative abundance in the sample > 0", file = file.path(resultsPath,vandamtaxafile), append = TRUE)

  }
  # begin loop over different variables

  ls_labels = list(
    salinity = c(
      "VD Salinity 1",
      "VD Salinity 2",
      "VD Salinity 3",
      "VD Salinity 4",
      "VD Salinity Indet",
      "VD Salinity Taxa used"
    ),
    heterotrophy = c(
      "VD N-Het 1",
      "VD N-Het 2",
      "VD N-Het 3",
      "VD N-Het 4",
      "VD N-Het Indet",
      "VD N-Het Taxa used"
    ),
    oxygen = c(
      "VD Oxygen 1",
      "VD Oxygen 2",
      "VD Oxygen 3",
      "VD Oxygen 4",
      "VD Oxygen 5",
      "VD Oxygen Indet",
      "VD Oxygen Taxa used"
    ),
    saprobity = c(
      "VD Saprobity 1",
      "VD Saprobity 2",
      "VD Saprobity 3",
      "VD Saprobity 4",
      "VD Saprobity 5",
      "VD Saprobity Indet",
      "VD Saprobity Taxa used"
    ),
    aero = c(
      "VD Aero 1",
      "VD Aero 2",
      "VD Aero 3",
      "VD Aero 4",
      "VD Aero 5",
      "VD Aero Indet",
      "VD Aero Taxa used"
    ),
    trophic = c(
      "VD Trophic 1",
      "VD Trophic 2",
      "VD Trophic 3",
      "VD Trophic 4",
      "VD Trophic 5",
      "VD Trophic 6",
      "VD Trophic 7",
      "VD Trophic Indet",
      "VD Trophic Taxa used"
    )
  )
  # list to hold loop outputs
  ls_output = list()

  for (i in 1:6){

    # setting up switches
    vdam_var = switch (i,
                       taxaInRA[, "vdams"],
                       taxaInRA[, "vdamnh"],
                       taxaInRA[, "vdamo2"],
                       taxaInRA[, "vdamsap"],
                       taxaInRA[, "vdamaero"],
                       taxaInRA[, "vdamtrop"])
    mess_var = switch(i,
                      "salinity",
                      "N-heterotrophy",
                      "oxygen requirements",
                      "saprobity",
                      "aero",
                      "trophic state")

    print(paste("Calculating Van Dam", mess_var))

    vdam_var[is.na(vdam_var)] = 0

    lp_data =  suppressWarnings(data.table(
      vd1 = unlist(taxaInRA[which(vdam_var == 1),
                            lapply(.SD, sum, na.rm = TRUE),
                            .SDcols = 1:(lastcol - 1)]),
      vd2 = unlist(taxaInRA[which(vdam_var == 2),
                            lapply(.SD, sum, na.rm = TRUE),
                            .SDcols = 1:(lastcol - 1)]),
      vd3 = unlist(taxaInRA[which(vdam_var == 3),
                            lapply(.SD, sum, na.rm = TRUE),
                            .SDcols = 1:(lastcol - 1)]),
      vd4 = unlist(taxaInRA[which(vdam_var == 4),
                            lapply(.SD, sum, na.rm = TRUE),
                            .SDcols = 1:(lastcol - 1)])
    ))
    if (i > 2){
      lp_data[, vd5 := unlist(taxaInRA[which(vdam_var == 5),
                                       lapply(.SD, sum, na.rm = TRUE),
                                       .SDcols = 1:(lastcol - 1)])]
    }

    if (i == 6){
      lp_data[, vd6 := unlist(taxaInRA[which(vdam_var == 6),
                                       lapply(.SD, sum, na.rm = TRUE),
                                       .SDcols = 1:(lastcol - 1)])]
      lp_data[, vd7 := unlist(taxaInRA[which(vdam_var == 7),
                                       lapply(.SD, sum, na.rm = TRUE),
                                       .SDcols = 1:(lastcol - 1)])]
    }

   #remove possible NAs
    lp_data[is.na(lp_data)] <- 0

    ## indet ##
    if (i < 3) {
      lp_data[      , v1 := round(100 - (vd1 + vd2 + vd3 + vd4), 1)]
    } else if (i == 4 | i == 5){
      lp_data[      , v1 := round(100 - (vd1 + vd2 + vd3 + vd4 + vd5), 1)]
    } else if (i == 6){
      lp_data[      , v1 := round(100 - (vd1 + vd2 + vd3 + vd4 + vd5 + vd6 +vd7), 1)]
    }
    lp_data[v1 < 0, v1 := 0]
    ##  taxa used ##
    # - prepare numeric column
    lp_data[      , v2 := numeric(lastcol - 1)]
    # - count the number of entries for each column that are not zero and
    # - are also not zero included in the vector of classes (i.e. vdams_v).
    te2 = purrr::map(.x = 1:(lastcol - 1),
                     .f = ~ taxaInRa_samples[, .x, with = F] > 0 &
                       vdam_var != 0)
    # - te2 is a list with one one dimensional array per element.
    # - each array represents one columns of taxaInRa
    # - We want to unlist the list and store it in a data.frame
    te3 = matrix(unlist(te2), ncol = lastcol - 1, byrow = TRUE) %>%
      as.data.frame()
    # rename columns
    names(te3) = names(taxaInRa_samples)
    # - until now the column of te3 are boolean - sum up to obtain number of
    # - species
    lp_data[, v2 := colSums(te3)]
    names(lp_data) = ls_labels[[i]]
    ls_output[[i]] = lp_data

    # write used taxa to file

    if (vandamReports & exists("resultsPath")){
      print_taxa = lapply(te2, function(x)taxaInRA$recognizedSp[x])
      for (print.i in seq_along(print_taxa)){
        write(paste(str_to_upper(mess_var), "- Sample:", colnames(taxaInRA)[print.i]),
              file.path(resultsPath, vandamtaxafile),
              append = TRUE)
        write(print_taxa[[print.i]], file.path(resultsPath, vandamtaxafile), append = TRUE)
      }
    }
    rm(vdam_var, te3, te2, lp_data, mess_var)
  }


  if (nrow(ls_output[[1]]) > 0){
    vandam.results <- data.frame(ls_output[[1]])
  }
  if (nrow(ls_output[[2]]) > 0){
    vandam.results <- cbind(vandam.results, data.frame(ls_output[[2]]))
  }
  if (nrow(ls_output[[3]]) > 0){
    vandam.results <- cbind(vandam.results, data.frame(ls_output[[3]]))
  }
  if (nrow(ls_output[[4]]) > 0){
    vandam.results <- cbind(vandam.results, data.frame(ls_output[[4]]))
  }
  if (nrow(ls_output[[5]]) > 0){
    vandam.results <- cbind(vandam.results, data.frame(ls_output[[5]]))
  }
  if (nrow(ls_output[[6]]) > 0){
    vandam.results <- cbind(vandam.results, data.frame(ls_output[[6]]))
  }

  vandam.results[is.na(vandam.results)] <- 0 #have to remove NAs for plotting
  return(vandam.results)


}
