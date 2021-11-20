#' Loads the 'Diat.Barcode' database into DiaThor in the correct format
#' @description
#' The package downloads and installs a wrapper for the 'Diat.Barcode' project. Besides citing the DiaThor package, the Diat.Barcode project should also be cited, as follows:
#' \itemize{
#' \item Rimet F., Gusev E., Kahlert M., Kelly M., Kulikovskiy M., Maltsev Y., Mann D., Pfannkuchen M., Trobajo R., Vasselon V., Zimmermann J., Bouchez A., 2019. Diat.barcode, an open-access curated barcode library for diatoms. Scientific Reports. https://www.nature.com/articles/s41598-019-51500-6
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @importFrom stringdist stringdist ain
#' @export diat_getDiatBarcode

diat_getDiatBarcode <- function() {

  ########## LINK WITH DIAT.BARCODE DATABASE
  #internal 'Diat.Barcode' number, when no internet connection is available
  intversion <- 10.1
  #get version number of latest 'Diat.Barcode'
  dic <- read.csv("http://www.francoiskeck.fr/work/diatbarcode/dic_version.csv", header = TRUE, stringsAsFactors = FALSE)
  if (exists("dic")){
    #is able to check the version
    version <- dic[dic$Flavor =="original",]
    version <- version$Version[which.max(as.numeric(as.POSIXlt(version$Date, format = "%d-%m-%Y")))]

    #compare both. If updates are needed, attempt them
    if (version == intversion){
      #updates are not needed
      print ("No updates needed for the 'Diat.barcode' database. Proceeding!")
      dbc <- diathor::dbc_offline
    } else {
      #updates are needed
      ########--------  Diat.Barcode download attempt. If it fails, tries to use internal database
      ## WARNING: CRAN package does not auto-update the Diat.Barcode database
      print("The diatom database in DiaThor is out of date")

      ###### THIS SECTION IS FOR THE CRAN PROJECT ONLY
      print("The CRAN version of the package does not auto-update the internal database. But the GitHub version does!")
      print("Using internal database, 'Diat.barcode' v.10.1 published on 25-06-2021")
      dbc <- diathor::dbc_offline
      ###### END OF CRAN VERSION

      ###### THIS SECTION IS FOR THE GITHUB PROJECT ONLY
      #  print("Attempting to download diat.barcode from website")
      # #
      #  #try to get the updated package
      #  tryCatch(
      #    {
      #     dbc <- diatbarcode::get_diatbarcode(version = "last") #loads the latest version of diat.barcode
      #      print("Latest version of Diat.barcode succesfully downloaded. Remember to credit accordingly!")
      #    },
      #    error = function(cond){
      #      print("Latest version of Diat.barcode cannot be downloaded")
      #      print("Using internal database, Diat.barcode v.10.1 published on 25-06-2021. It might need to be updated")
      #      dbc <- diathor::dbc_offline
      #
      #    }
      #  )
      ###### END OF GITHUB VERSION

    }
  } else {
    print("Latest version of 'Diat.barcode' unknown")
    print("Using internal database, 'Diat.barcode' v.10.1 published on 25-06-2021")
    dbc <- diathor::dbc_offline
  }
  ### Double checks that database got loaded correctly or cancels alltogether
  if (exists("dbc")){
    #database loaded ok

  } else { #if everything fails, cancels
    stop("Database could not be downloaded or retrieved from internal package. Cancelling calculations")
  }
  ########## END LINK WITH DIAT.BARCODE DATABASE
  #Remove duplicate by field "species" in diat.barcode
  dbc2 <- as.data.frame(dbc[!duplicated(dbc[,"species"]),]) #transforms dbc to a dataframe
  ecodata <- dbc2[which(colnames(dbc2)=="species" ):ncol(dbc2)] #keeps only from the "species" column onwards, to keep the ecological data
  #Update the internal taxa list
  taxaList <- diathor::taxaList
  inner_taxaList <- diat_taxaList()
  taxaList <- as.data.frame(c(inner_taxaList$species,
                              ecodata$species))
  colnames(taxaList) <- "species"

  #Returns ecodata and taxon list
  return (list(ecodata, taxaList))


}
