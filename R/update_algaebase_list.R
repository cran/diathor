#' Updates a list of taxa names against AlgaeBase.org
#' @description
#' The package is a wrapper for the algaeClassify package from the 'AlgaeBase.org' database. Besides citing the DiaThor package, the Diat.Barcode project should also be cited, as follows:
#' \itemize{
#' \item Patil, V., Seltmann, T., Salmaso, N., Anneville, O., Lajeunesse, M., Straile, D. 2019. algaeClassify: Determine Phytoplankton Functional Groups Based on Functional Traits.
#' }
#' @param listToUpdate the name of the taxa (genus, species, variety) to be checked against the internal DB
#' @param overwrite if overwrite = T (default), if a new name is found, it replaces the original one. The corrected list is returned. Otherwise, a new list is returned, with the original taxon and the new name found
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @import algaeClassify
#' @export update_algaebase_list

#UPDATE NAMES FROM ALGAEBASE
update_algaebase_list <- function(listToUpdate, overwrite = T){
  if (overwrite == T) {
    pb <- txtProgressBar(min = 1, max = nrow(listToUpdate), style = 3)
    for (i in 1:nrow(listToUpdate)) {
      oriRow <- listToUpdate[i,]
      oldname <- oriRow[,"fullspecies"]
      first_word <- strsplit(oldname, " ")[[1]][1] #first word
      other_words <- strsplit(oldname, " ")[[1]][-1] #other words
      other_words <- paste(other_words, collapse = " ") #other words
      #search name in algaebase
      new_name_res <- algae_search(genus=first_word, species=other_words, long=FALSE)
      #if it finds a match
      if (new_name_res$exact.match[1] == 1){
        new_genus <- new_name_res$genus
        new_species <- new_name_res$species
        #creates new species
        new_name <- paste(new_genus, new_species, sep = " ")
        oriRow[,"fullspecies"] <- new_name
        listToUpdate[i,] <- oriRow
      }
      setTxtProgressBar(pb, i)
    }
    #close progressbar
    close(pb)

    #trimws
    listToUpdate$fullspecies <- trimws(listToUpdate$fullspecies)
    #remove duplicate species
    listToUpdate <- as.data.frame(listToUpdate[!duplicated(listToUpdate[,"fullspecies"]),])

    return(listToUpdate)

  } else {
    # no overwrite, new column
    listToUpdate$new_name <- NA
    pb <- txtProgressBar(min = 1, max = nrow(listToUpdate), style = 3)
    for (i in 1:nrow(listToUpdate)) {
      oriRow <- listToUpdate[i,]
      oldname <- oriRow[,"fullspecies"]
      first_word <- strsplit(oldname, " ")[[1]][1] #first word
      other_words <- strsplit(oldname, " ")[[1]][-1] #other words
      other_words <- paste(other_words, collapse = " ") #other words
      #search name in algaebase
      new_name_res <- algae_search(genus=first_word, species=other_words, long=FALSE)
      #if it finds a match
      if (new_name_res$exact.match[1] == 1){
        new_genus <- new_name_res$genus
        new_species <- new_name_res$species
        #creates new species
        new_name <- paste(new_genus, new_species, sep = " ")
        oriRow[,"new_name"] <- new_name
        listToUpdate[i,] <- oriRow
      }
      setTxtProgressBar(pb, i)
    }
    #close progressbar
    close(pb)

    #trimws
    listToUpdate$fullspecies <- trimws(listToUpdate$fullspecies)
    #remove duplicate species
    listToUpdate <- as.data.frame(listToUpdate[!duplicated(listToUpdate[,"fullspecies"]),])

    return(listToUpdate)
  }

}
