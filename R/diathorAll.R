#' Runs all the DiaThor functions in a pipeline
#' @param species_df The data frame with your species data. Species as rows, Sites as columns. If empty, a dialog will prompt for a CSV file
#' @param isRelAb Boolean. If set to 'TRUE' it means that your species' data is the relative abundance of each species per site. If FALSE, it means that it the data corresponds to absolute densities. Default = FALSE
#' @param maxDistTaxa Integer. Number of characters that can differ in the species' names when compared to the internal database's name in the heuristic search. Default = 2
#' @param resultsPath String. Path to the output folder. If none specified (default), a dialog box will prompt to select it
#' @param calculateguilds Boolean. If set to 'TRUE' the percentage of abundance of each diatom guild will be calculated. Default = TRUE
#' @param vandam Boolean. If set to 'TRUE' the Van Dam classifications will be calculated in the Output. Default = TRUE
#' @param vandamReports Boolean. If set to 'TRUE' the detailed reports for the Van Dam classifications will be reported in the Output. Default = TRUE
#' @param singleResult Boolean. If set to 'TRUE' all results will go into a single output file. If FALSE, separate output files will be generated. Default = TRUE
#' @param plotAll Boolean. If set to 'TRUE', plots will be generated for each Output in a PDF file. Default = TRUE
#' @param exportFormat Integer. If = 1: only a CSV (external file) will be generated with the output matrices; 2: only an internal R dataframe will be generated; 3: both a CSV and an internal R dataframe are generated. Default = 3
#' @param exportName String. Prefix for the CSV exported file. Default = "DiaThor_results"
#' @param color Color code (hex). Default color for bar charts and lolipop plots. Default = "#0073C2"
#' @description
#' The diaThorAll function is the master function of the package. It calculates all outputs from the data, and places them in the Output folder
#' The input file for the package is a dataframe or an external CSV file. Species should be listed as rows, with species' names in column 1 (column name should be "species")
#' The other columns (samples) have to contain the abundance of each species (relative or absolute) in each sample.
#' The first row of the file has to contain the headers with the sample names. Remember that a column named "species" is mandatory, containing the species' names
#' If a dataframe is not specified as a parameter (species_df), the package will show a dialog box to search for the CSV file
#' A second dialog box will help set up an Output folder, where all outputs from the package will be exported to (dataframes, CSV files, plots in PDF)
#' The package also downloads and installs a wrapper for the 'Diat.Barcode' project. Besides citing the DiaThor package, the Diat.Barcode project should also be cited, as follows:
#' \itemize{
#' \item Rimet, Frederic; Gusev, Evgenuy; Kahlert, Maria; Kelly, Martyn; Kulikovskiy, Maxim; Maltsev, Yevhen; Mann, David; Pfannkuchen, Martin; Trobajo, Rosa; Vasselon, Valentin; Zimmermann, Jonas; Bouchez, Agnès. 2018. "Diat.barcode, an open-access barcode library for diatoms". Scientific Reports,9, 15116. https://doi.org/10.15454/TOMBYZ
#' }
#' Sample data in the examples is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951. https://doi.org/10.1016/j.ecolind.2019.105951
#' }
#' @examples
#' \donttest{
#' # Example using sample data included in the package (sampleData):
#' data("diat_sampleData")
#' # In the example, a temporary directory will be used in resultsPath
#' allResults <- diaThorAll(diat_sampleData, resultsPath = tempdir())
#' }
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @importFrom grDevices dev.off pdf
#' @importFrom tidyr gather
#' @import ggplot2
#' @import utils
#' @import data.table
#' @rawNamespace import(purrr, except = transpose)
#' @import stringr
#' @importFrom tibble tibble
#' @export diaThorAll

###### ---------- MASTER FUNCTION, CALCULATES EVERYTHING WITH ALL POSSIBLE OUTPUTS BY DEFAULT  ---------- ########

diaThorAll <- function(species_df, isRelAb=FALSE, maxDistTaxa = 2, resultsPath, calculateguilds = TRUE, vandam = TRUE, vandamReports = TRUE, singleResult = TRUE, exportFormat = 3, exportName = "DiaThor_results", plotAll = TRUE, color = "#0073C2"){
  resultmat <- diat_loadData(species_df, isRelAb, maxDistTaxa, resultsPath)
  morpho.results <- diat_morpho(resultmat, isRelAb)

  if (exists("morpho.results")){
    numcloroplastos.result <- morpho.results[[1]]
    shpcloroplastos.result <- morpho.results[[2]]
    biovol.val.result <- morpho.results[[3]]
  }
  size.results <- diat_size(resultmat)
  guilds.results <- diat_guilds(resultmat)
  diversity.results <- diat_diversity(resultmat)
  vandam.results <- diat_vandam(resultmat)
  ips.results <- diat_ips(resultmat)
  tdi.results <- diat_tdi(resultmat)
  idp.results <- diat_idp(resultmat)
  ilm.results <- diat_ilm(resultmat)
  des.results <- diat_des(resultmat)
  epid.results <- diat_epid(resultmat)
  idap.results <- diat_idap(resultmat)
  idch.results <- diat_idch(resultmat)
  lobo.results <- diat_lobo(resultmat)
  sla.results <- diat_sla(resultmat)
  spear.results <- diat_spear(resultmat)

  #sampledata
  sampleNames <- resultmat[[3]]
  resultsPath <- resultmat[[4]]

  ########### RESULTS TABLES ############
  if (singleResult == TRUE){ #by default exports a single spreadsheet with all the results
    singleTable <- NULL
    singleTable <- data.frame(c(if(exists("diversity.results")){diversity.results},
                     if(exists("numcloroplastos.result")){as.data.frame(numcloroplastos.result)},
                     if(exists("shpcloroplastos.result")){shpcloroplastos.result},
                     if(exists("biovol.val.result")){biovol.val.result},
                     if(exists("size.results")){size.results},
                     if(exists("guilds.results")){guilds.results},
                     if(exists("vandam.results")){vandam.results},
                     if(exists("ips.results")){ips.results},
                     if(exists("tdi.results")){tdi.results},
                     if(exists("idp.results")){idp.results},
                     if(exists("ilm.results")){ilm.results},
                     if(exists("des.results")){des.results},
                     if(exists("epid.results")){epid.results},
                     if(exists("idap.results")){idap.results},
                     if(exists("idch.results")){idch.results},
                     if(exists("lobo.results")){lobo.results},
                     if(exists("sla.results")){sla.results},
                     if(exists("spear.results")){spear.results}
                     ))
    #removes the Precision columns
    #singleTable[ , -which(names(singleTable) %in% "Precision")]
    singleTable <- singleTable[ , -which(startsWith(names(singleTable),"Precision")) ]

    rownames(singleTable) <- sampleNames
  } else { #separate files for each result
    listOfTables <- list(if(exists("diversity.results")){diversity.results},
                     if(exists("numcloroplastos.result")){numcloroplastos.result},
                     if(exists("shpcloroplastos.result")){shpcloroplastos.result},
                     if(exists("biovol.val.result")){biovol.val.result},
                     if(exists("size.results")){size.results},
                     if(exists("guilds.results")){guilds.results},
                     if(exists("vandam.results")){vandam.results},
                     if(exists("ips.results")){ips.results},
                     if(exists("tdi.results")){tdi.results},
                     if(exists("idp.results")){idp.results},
                     if(exists("ilm.results")){ilm.results},
                     if(exists("des.results")){des.results},
                     if(exists("epid.results")){epid.results},
                     if(exists("idap.results")){idap.results},
                     if(exists("idch.results")){idch.results},
                     if(exists("lobo.results")){lobo.results},
                     if(exists("sla.results")){sla.results},
                     if(exists("spear.results")){spear.results}
    )

    names(listOfTables) <- c(if(exists("diversity.results")){"Diversity"},
                         if(exists("numcloroplastos.result")){"Chloroplast number"},
                         if(exists("shpcloroplastos.result")){"Chloroplast shape"},
                         if(exists("biovol.val.result")){"Biovolume"},
                         if(exists("size.results")){"Size classes"},
                         if(exists("guilds.results")){"Guilds"},
                         if(exists("vandam.results")){"VanDam"},
                         if(exists("ips.results")){"IPS index"},
                         if(exists("tdi.results")){"TDI index"},
                         if(exists("idp.results")){"IDP index"},
                         if(exists("ilm.results")){"ILM index"},
                         if(exists("des.results")){"DES index"},
                         if(exists("epid.results")){"EPID index"},
                         if(exists("idap.results")){"IDAP index"},
                         if(exists("idch.results")){"IDCH index"},
                         if(exists("lobo.results")){"LOBO index"},
                         if(exists("sla.results")){"SLA index"},
                         if(exists("spear.results")){"SPEAR index"}
                         )

  }


  ########### END RESULTS TABLES

  ########### START PLOTS ############

  loli.plot <- function(result, ylabel, ylow, yhigh, samplenames, color = "#0073C2"){
    x <- rownames(result)
    data <- result
    data[is.na(data)] = 0

    for(i in 1:ncol(result)) {
      y <- data[,i]

      return(ggplot2::ggplot(data, aes(x=samplenames, y=y)) +
              geom_segment( aes(x=samplenames, xend=samplenames, y=0, yend=y), color="grey") +
              geom_point( color=color, size=4) +
              theme_light() +
              theme(
                panel.grid.major.x = element_blank(),
                panel.border = element_blank(),
                axis.ticks.x = element_blank()
              ) +
              ylim(ylow, yhigh) +
              xlab("") +
              ylab(ylabel)
      )
    }
  }

  percentbarchart.plot <- function(result, title){
    x <- colnames(result) #checks if taxa used colum exists and removes it
    x <- substr(x, nchar(x)-8, nchar(x)) == 'Taxa used' #checks if taxa used colum exists and removes it
    if (tail(x, n=1)==TRUE) {
      result <- result[,-ncol(result)] #and removes it
    }
    sampleCol <- rep(sampleNames, ncol(result)) #gets sample names
    result <- tidyr::gather(result) #uses tidyr to rearrange the dataframe in a single column
    result$sampleCol <- sampleCol #adds another column with the sample names
    colors <- c("#CC1C00", "#5C88DA", "#84BD00", "#FFCD00", "#7C878E", "#E64B35", "#4DBBD5", "#01A087", "#3C5488", "#F39B7F", "#FF410D99", "#6EE2FF99", "#F7C53099", "#95CC5E99", "#D0DFE699", "#F79D1E99", "#748AA699", "#82451c", "#4b7751", "#5fa413" )
    key <- result$key
    value <- result$value
    print(ggplot2::ggplot(result, aes(fill=key, y=value, x=sampleCol)) +
            geom_bar(position="fill", stat="identity") +
            scale_fill_manual(values=colors) +
            xlab("Samples") +
            ylab(title)) #graphs
  }

  result.plots.allToPDF <- function(color ="#0073C2"){
    print("Exporting all plots to PDF, please wait...")
    # Open a pdf file
    pdf(file.path(resultsPath, "Plots.pdf"))
#   pdf(paste(resultsPath, "\\", "Plots2.pdf", sep=""))

    #Plots all resulting graphs (if exist)
    if(exists("diversity.results")){

      print(loli.plot(as.data.frame(diversity.results[,1]), "Species richess", 0, max(as.data.frame(diversity.results[,1])), samplenames=rownames(diversity.results)))
      print(loli.plot(as.data.frame(diversity.results[,2]), "Shannon's diversity", 0, max(as.data.frame(diversity.results[,2])), samplenames=rownames(diversity.results)))
      print(loli.plot(as.data.frame(diversity.results[,3]), "Evenness", 0, max(as.data.frame(diversity.results[,3])), samplenames=rownames(diversity.results)))
    }
    if(exists("biovol.val.result")){
      print(loli.plot(as.data.frame(biovol.val.result[,1]), "Biovolume", 0, max(as.data.frame(biovol.val.result[,1])), samplenames=rownames(biovol.val.result)))
    }
    if(exists("numcloroplastos.result")){
      percentbarchart.plot(numcloroplastos.result, "Number of chloroplasts") #default: piled bars
    }
    if(exists("shpcloroplastos.result")){
      percentbarchart.plot(shpcloroplastos.result, "Shape of chloroplasts") #default: piled bars
    }
    if(exists("size.results")){
      percentbarchart.plot(size.results, "Size classes") #default: piled bars
    }
      if(exists("guilds.results")){
      percentbarchart.plot(guilds.results, "Guilds")#default: piled bars
    }
    if(exists("vandam.results")){
      vdamSalinity <- vandam.results[,startsWith(colnames(vandam.results),"VD.Salinity")]
      vdamNHeterotrophy <- vandam.results[,startsWith(colnames(vandam.results),"VD.N.Het")]
      vdamOxygen <- vandam.results[,startsWith(colnames(vandam.results),"VD.Oxygen")]
      vdamSaprobity <- vandam.results[,startsWith(colnames(vandam.results),"VD.Saprobity")]
      vdamAero <- vandam.results[,startsWith(colnames(vandam.results),"VD.Aero")]
      vdamTrophic <- vandam.results[,startsWith(colnames(vandam.results),"VD.Trophic")]
      #remove the Taxa Used column
      if (ncol(vdamSalinity) > 0){
        vdamSalinity <- vdamSalinity[1:(ncol(vdamSalinity)-1)]
      }
      if (ncol(vdamNHeterotrophy) > 0){
        vdamNHeterotrophy <- vdamNHeterotrophy[1:(ncol(vdamNHeterotrophy)-1)]
      }
      if (ncol(vdamOxygen) > 0){
        vdamOxygen <- vdamOxygen[1:(ncol(vdamOxygen)-1)]
      }
      if (ncol(vdamSaprobity) > 0){
        vdamSaprobity <- vdamSaprobity[1:(ncol(vdamSaprobity)-1)]
      }
      if (ncol(vdamAero) > 0){
        vdamAero <- vdamAero[1:(ncol(vdamAero)-1)]
      }
      if (ncol(vdamTrophic) > 0){
        vdamTrophic <- vdamTrophic[1:(ncol(vdamTrophic)-1)]
      }
    }
    if(exists("vdamSalinity")){
      if (ncol(vdamSalinity)>0){
        percentbarchart.plot(vdamSalinity, "Salinity")
      }
    }
    if(exists("vdamNHeterotrophy")){
      if (ncol(vdamNHeterotrophy)>0){
        percentbarchart.plot(vdamNHeterotrophy, "N-Heterotrophy")
      }
    }
    if(exists("vdamOxygen")){
      if (ncol(vdamOxygen)>0){
        percentbarchart.plot(vdamOxygen, "Oxygen preferences")
      }
    }
    if(exists("vdamSaprobity")){
      if (ncol(vdamSaprobity)>0){
        percentbarchart.plot(vdamSaprobity, "Saprobity")
      }
    }
    if(exists("vdamAero")){
      if (ncol(vdamAero)>0){
        percentbarchart.plot(vdamAero, "Moisture")
      }
    }
    if(exists("vdamTrophic")){
      if (ncol(vdamTrophic)>0){
        percentbarchart.plot(vdamTrophic, "Trophic state")
      }
    }
    if(exists("ips.results")){
      print(loli.plot(as.data.frame(ips.results[,1]), "IPS", 0, 5, samplenames=rownames(ips.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #raw index. 0 = NA
      print(loli.plot(as.data.frame(ips.results[,2]), "IPS - Standardized", 0, 20, samplenames=rownames(ips.results))) #standard 20
    }
    if(exists("tdi.results")){
      print(loli.plot(as.data.frame(tdi.results[,1]), "TDI - Standardized", 0, 20, samplenames=rownames(tdi.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #standard 20. 0 = NA
      print(loli.plot(as.data.frame(tdi.results[,2]), "TDI - %", 0, 100, samplenames=rownames(tdi.results))) #%
    }
    if(exists("idp.results")){
      print(loli.plot(as.data.frame(idp.results[,1]), "IDP", 0, 4, samplenames=rownames(idp.results))) #raw index
      print(loli.plot(as.data.frame(idp.results[,2]), "IDP - Standardized", 0, 20, samplenames=rownames(idp.results))) #standard 20
    }
    if(exists("ilm.results")){
      print(loli.plot(as.data.frame(ilm.results[,1]), "ILM", 0, 5, samplenames=rownames(ilm.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #raw index. 0 = NA
      print(loli.plot(as.data.frame(ilm.results[,2]), "ILM - Standardized", 0, 20, samplenames=rownames(ilm.results))) #standard 20
    }
    if(exists("des.results")){
      print(loli.plot(as.data.frame(des.results[,1]), "DES", 0, 5, samplenames=rownames(des.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #raw index. 0 = NA
      print(loli.plot(as.data.frame(des.results[,2]), "DES - Standardized", 0, 20, samplenames=rownames(des.results))) #standard 20
    }
    if(exists("epid.results")){
      print(loli.plot(as.data.frame(epid.results[,1]), "EPID", 0, 4, samplenames=rownames(epid.results))) #raw index
      print(loli.plot(as.data.frame(epid.results[,2]), "EPID - Standardized", 0, 20, samplenames=rownames(epid.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #standard 20. 0 = NA
    }
    if(exists("idap.results")){
      print(loli.plot(as.data.frame(idap.results[,1]), "IDAP", 0, 5, samplenames=rownames(idap.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #raw index. 0 = NA
      print(loli.plot(as.data.frame(idap.results[,2]), "IDAP - Standardized", 0, 20, samplenames=rownames(idap.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #standard 20. 0 = NA
    }
    if(exists("idch.results")){
      print(loli.plot(as.data.frame(idch.results[,1]), "IDCH", 0, 8, samplenames=rownames(idch.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #raw index. 0 = NA
      print(loli.plot(as.data.frame(idch.results[,2]), "IDCH - Standardized", 0, 20, samplenames=rownames(idch.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #standard 20. 0 = NA
    }
    if(exists("lobo.results")){
      print(loli.plot(as.data.frame(lobo.results[,1]), "LOBO", 0, 4, samplenames=rownames(lobo.results))+ geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #raw index. 0 = NA
      print(loli.plot(as.data.frame(lobo.results[,2]), "LOBO - Standardized", 0, 20, samplenames=rownames(lobo.results))) #standard 20
    }
    if(exists("sla.results")){
      print(loli.plot(as.data.frame(sla.results[,1]), "SLA", 0, 4, samplenames=rownames(sla.results))) #raw index
      print(loli.plot(as.data.frame(sla.results[,2]), "SLA - Standardized", 0, 20, samplenames=rownames(sla.results)) + geom_hline(yintercept=1, linetype="dashed", color = "darkgray", size=1)) #standard 20. 0 = NA
    }
    if(exists("spear.results")){
      print(loli.plot(as.data.frame(spear.results[,1]), "SPEAR", 0, 100, samplenames=rownames(spear.results))) #raw index
    }

    # Close the pdf file
    dev.off()
    print("Plots exported!")

  }

  if (plotAll == TRUE){
    result.plots.allToPDF()
  }

  ########### END PLOTS

  #EXPORT AS CSV
  if (exportFormat == 1) {
    if (singleResult == TRUE) {
      filename <- paste(exportName, " - Results.csv", sep ="")
      write.csv(singleTable, file.path(resultsPath, filename))

    } else {
      for (i in seq_along(listOfTables)) {
        filename <- paste(exportName, " - ",names(listOfTables)[i], ".csv")
        write.csv(singleTable, file.path(resultsPath, filename))
      }
    }
  }

  #EXPORT AS INTERNAL DATAFRAME
  if (exportFormat == 2) {

    if (singleResult == TRUE) {
      return(singleTable)
    } else {
      resultList <- list(as.data.frame(listOfTables[[1]]),
                         as.data.frame(listOfTables[[2]]),
                         as.data.frame(listOfTables[[3]]),
                         as.data.frame(listOfTables[[4]]),
                         as.data.frame(listOfTables[[5]]),
                         as.data.frame(listOfTables[[6]]),
                         as.data.frame(listOfTables[[7]]),
                         as.data.frame(listOfTables[[8]]),
                         as.data.frame(listOfTables[[9]]),
                         as.data.frame(listOfTables[[10]]),
                         as.data.frame(listOfTables[[11]]),
                         as.data.frame(listOfTables[[12]]),
                         as.data.frame(listOfTables[[13]]),
                         as.data.frame(listOfTables[[14]]),
                         as.data.frame(listOfTables[[15]]),
                         as.data.frame(listOfTables[[16]]),
                         as.data.frame(listOfTables[[17]]),
                         as.data.frame(listOfTables[[18]])
                          )
      names(resultList) <- c("Diversity",
                             "ChloroplastNumber",
                             "ChloroplastShape",
                             "Biovolume",
                             "SizeClasses",
                             "Guilds",
                             "VanDam",
                             "IPS",
                             "TDI",
                             "IDP",
                             "ILM",
                             "DES",
                             "EPID",
                             "EPID",
                             "IDAP",
                             "IDCH",
                             "LOBO",
                             "SPEAR"
                             )
      return(resultList)
    }
  }

  #EXPORT AS BOTH CSV AND INTERNAL DATAFRAME - Default
  if (exportFormat == 3) {
    if (singleResult == TRUE) {
      filename <- paste(exportName, " - Results.csv", sep ="")
      # write.csv(singleTable, paste(resultsPath, "\\", filename, sep=""))
      write.csv(singleTable, file.path(resultsPath, filename))

      return(singleTable)
    } else {
      for (i in seq_along(listOfTables)) {
        filename <- paste(exportName, " - ",names(listOfTables)[i], ".csv", sep ="")
        # write.csv(listOfTables[[i]], paste(resultsPath, "\\", filename, sep=""))
        write.csv(singleTable, file.path(resultsPath, filename))

      }

      resultList <- list(as.data.frame(listOfTables[[1]]),
           as.data.frame(listOfTables[[2]]),
           as.data.frame(listOfTables[[3]]),
           as.data.frame(listOfTables[[4]]),
           as.data.frame(listOfTables[[5]]),
           as.data.frame(listOfTables[[6]]),
           as.data.frame(listOfTables[[7]]),
           as.data.frame(listOfTables[[8]]),
           as.data.frame(listOfTables[[9]]),
           as.data.frame(listOfTables[[10]]),
           as.data.frame(listOfTables[[11]]),
           as.data.frame(listOfTables[[12]]),
           as.data.frame(listOfTables[[13]]),
           as.data.frame(listOfTables[[14]]),
           as.data.frame(listOfTables[[15]]),
           as.data.frame(listOfTables[[16]]),
           as.data.frame(listOfTables[[17]]),
           as.data.frame(listOfTables[[18]])
      )
      names(resultList) <- c("Diversity",
                             "ChloroplastNumber",
                             "ChloroplastShape",
                             "Biovolume",
                             "SizeClasses",
                             "Guilds",
                             "VanDam",
                             "IPS",
                             "TDI",
                             "IDP",
                             "ILM",
                             "DES",
                             "EPID",
                             "EPID",
                             "IDAP",
                             "IDCH",
                             "LOBO",
                             "SPEAR"
      )
      return(resultList)
    }
  }


}
