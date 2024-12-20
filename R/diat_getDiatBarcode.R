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
  # Internal 'Diat.Barcode' version, when no internet connection is available
  intversion <- 10.1

  url <- "http://www.francoiskeck.fr/work/diatbarcode/dic_version.csv"

  # Attempt to download the latest version number of 'Diat.Barcode'
  # first we check if the file exists
  file_exists <- tryCatch(
    {
      # Try to open a connection to the URL
      suppressWarnings({
        con <- url(url, "r")
        close(con)


      })
      TRUE
    },
    error = function(e) {
      print("Error checking diatbarcode existence: 404, file not found")
      FALSE
    }
  )
  #if the file exists, we try to download it. This is extra steps are done to prevent 404 errors from aborting the whole function
  if(file_exists == T){
    dic <- tryCatch(
      {
        dic <- read.csv("http://www.francoiskeck.fr/work/diatbarcode/dic_version.csv", header = TRUE, stringsAsFactors = FALSE)
      },
      error = function(e) {
        print("Error occurred downloading diatbarcode")
      }
    )
  } else {

  }


  if (!exists("dic")) {
    # If dic is NULL, load internal database
    print("Latest version of Diat.barcode unknown.")
    print("Using internal database, 'Diat.barcode' v.10.1 published on 25-06-2021.")
    dbc <- diathor::dbc_offline
  } else {
    # If dic is not NULL, check the version
    version <- dic[dic$Flavor == "original",]
    version <- version$Version[which.max(as.numeric(as.POSIXlt(version$Date, format = "%d-%m-%Y")))]

    if (version == intversion) {
      # No updates needed
      print("No updates needed for the 'Diat.barcode' database. Proceeding!")
      dbc <- diathor::dbc_offline
    } else {
      # Updates are needed
      print("The diatom database in DiaThor is out of date.")

      ###### THIS SECTION IS FOR THE CRAN PROJECT ONLY
      # The CRAN version does not auto-update the internal database
      message("The CRAN version of the package does not auto-update the internal database. Using internal database, 'Diat.barcode' v.10.1 published on 25-06-2021.")
      dbc <- diathor::dbc_offline
      ###### END OF CRAN VERSION

      ###### THIS SECTION IS FOR THE GITHUB PROJECT ONLY
      # print("Attempting to download diat.barcode from website")
      # #
      # # try to get the updated package
      # dbc <- tryCatch(
      #   {
      #     diatbarcode::get_diatbarcode(version = "last") # loads the latest version of diat.barcode
      #   },
      #   error = function(e) {
      #     print("Latest version of Diat.barcode cannot be downloaded: ", e$message)
      #     print("Using internal database, 'Diat.barcode' v.10.1 published on 25-06-2021. It might need to be updated.")
      #     diathor::dbc_offline
      #   }
      #
      # )
      # updated <- isTRUE(!identical(dbc, diathor::dbc_offline))
      # if (updated==T){print("Latest version of Diat.barcode successfully downloaded. Remember to credit accordingly!")}

      ###### END OF GITHUB VERSION

    }
  }

  # Return ecodata and taxa list
  return(dbc)
}
