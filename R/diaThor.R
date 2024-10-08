#' DiaThor: A package to calculate multiple diatom-based biotic indices
#'
#' @description
#' The package calculates multiple biotic indices using diatoms from environmental samples. Diatom species are recognized by their species' name using a heuristic search, and their ecological data is retrieved from multiple sources.
#' Morphological information about the species is retrieved from the 'Diat.Barcode' project:
#' \itemize{
#' \item Rimet F., Gusev E., Kahlert M., Kelly M., Kulikovskiy M., Maltsev Y., Mann D., Pfannkuchen M., Trobajo R., Vasselon V., Zimmermann J., Bouchez A., 2019. Diat.barcode, an open-access curated barcode library for diatoms. Scientific Reports. https://www.nature.com/articles/s41598-019-51500-6
#' }
#' Size class classification is obtained from:
#' \itemize{
#' \item Rimet F. & Bouchez A., 2012. Life-forms, cell-sizes and ecological guilds of diatoms in European rivers. Knowledge and management of aquatic ecosystems, 406: 1-14. DOI:10.1051/kmae/2012018
#' }
#' Guild classification is obtained from:
#' \itemize{
#' \item Rimet F. & Bouchez A., 2012. Life-forms, cell-sizes and ecological guilds of diatoms in European rivers. Knowledge and management of aquatic ecosystems, 406: 1-14. DOI:10.1051/kmae/2012018
#' }
#' The combined classification of size classes and guilds is obtained from:
#' \itemize{
#' \item B-Béres, V., Török, P., Kókai, Z., Lukács, Á., Enikő, T., Tóthmérész, B., & Bácsi, I. (2017). Ecological background of diatom functional groups: Comparability of classification systems. Ecological Indicators, 82, 183-188.
#' }
#' Van Dam classification is obtained form:
#' \itemize{
#' \item Van Dam, H., Mertens, A., & Sinkeldam, J. (1994). A coded checklist and ecological indicator values of freshwater diatoms from the Netherlands. Netherland Journal of Aquatic Ecology, 28(1), 117-133.
#' }
#' Diversity index (Shannons H') is calculated using the vegan package, following:
#' \itemize{
#' \item Shannon, C. E., and Weaver, W. (1949). ‘The Mathematical Theory of Communication.’ (University of Illinios Press: Urbana, IL, USA.)
#' }
#' Species tolerance and their ecological information to calculate each biotic index is retrieved from their original sources:
#' \itemize{
#' \item IPS: Coste, M. (1982). Étude des méthodes biologiques d’appréciation quantitative de la qualité des eaux. Rapport Cemagref QE Lyon-AF Bassin Rhône Méditerranée Corse.
#' }
#' \itemize{
#' \item TDI: Kelly, M. G., & Whitton, B. A. (1995). The trophic diatom index: a new index for monitoring eutrophication in rivers. Journal of Applied Phycology, 7(4), 433-444.
#' }
#' \itemize{
#' \item IDP: Gómez, N., & Licursi, M. (2001). The Pampean Diatom Index (IDP) for assessment of rivers and streams in Argentina. Aquatic Ecology, 35(2), 173-181.
#' }
#' \itemize{
#' \item DES: Descy, J. P. 1979. A new approach to water quality estimation  using  diatom.  Beih.  Nov  Hedw. 64:305-323
#' }
#' \itemize{
#' \item EPID: Dell'Uomo, A. (1996). Assessment of water quality of an Apennine river as a pilot study for diatom-based monitoring of Italian watercourses. Use of algae for monitoring rivers, 65-72.
#' }
#' \itemize{
#' \item IDAP: Prygiel, J., & Coste, M. (1993). The assessment of water quality in the Artois-Picardie water basin (France) by the use of diatom indices. Hydrobiologia, 269(1), 343-349.
#' }
#' \itemize{
#' \item ID-CH: Hürlimann J., Niederhauser P. 2007: Méthodes d’analyse et d’appréciation des cours d’eau. Diatomées Niveau R (région). État de l’environnement n° 0740. Office fédéral de l’environnement, Berne. 132 p
#' }
#' \itemize{
#' \item ILM: Leclercq, L., & Maquet, B. (1987). Deux nouveaux indices diatomique et de qualité chimique des eaux courantes. Comparaison avec différents indices existants. Cahier de Biology Marine, 28, 303-310.
#' }
#' \itemize{
#' \item LOBO: Lobo, E. A., Callegaro, V. L. M., & Bender, E. P. (2002). Utilização de algas diatomáceas epilíticas como indicadoras da qualidade da água em rios e arroios da Região Hidrográfica do Guaíba, RS, Brasil. Edunisc.
#' }
#' \itemize{
#' \item LOBO: Lobo, E. A., Bes, D., Tudesque, L., & Ector, L. (2004). Water quality assessment of the Pardinho River, RS, Brazil, using epilithic diatom assemblages and faecal coliforms as biological indicators. Vie et Milieu, 54(2-3), 115-126.
#' }
#' \itemize{
#' \item SLA: Sládecek, V. (1986). Diatoms as indicators of organic pollution. Acta hydrochimica et hydrobiologica, 14(5), 555-566.
#' }
#' \itemize{
#' \item SPEAR(herbicides): Wood, R. J., Mitrovic, S. M., Lim, R. P., Warne, M. S. J., Dunlop, J., & Kefford, B. J. (2019). Benthic diatoms as indicators of herbicide toxicity in rivers–A new SPEcies At Risk (SPEARherbicides) index. Ecological Indicators, 99, 203-213.
#' }
#' \itemize{
#' \item PBIDW: Castro-Roa, D., & Pinilla-Agudelo, G. (2014). Periphytic diatom index for assessing the ecological quality of the Colombian Andean urban wetlands of Bogotá. Limnetica, 33(2), 297-312.
#' }
#' \itemize{
#' \item DISP: Stenger-Kovács, C., Körmendi, K., Lengyel, E., Abonyi, A., Hajnal, É., Szabó, B., Buczkó, K. & Padisák, J. (2018). Expanding the trait-based concept of benthic diatoms: Development of trait-and species-based indices for conductivity as the master variable of ecological status in continental saline lakes. Ecological Indicators, 95, 63-74.
#' }
#' Sample data included in the package is taken from:
#' \itemize{
#' \item Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge; Gómez, Nora. 2020. "Exploring the use of nuclear alterations, motility and ecological guilds in epipelic diatoms as biomonitoring tools for water quality improvement in urban impacted lowland streams". Ecological Indicators, 110, 105951.
#' }
#' @section Functions:
#' diat_loadData()
#' diat_morpho()
#' diat_size()
#' diat_diversity()
#' diat_guilds()
#' diat_vandam()
#' diat_loadData()
#' diat_ips()
#' diat_tdi()
#' diat_idp()
#' diat_des()
#' diat_epid()
#' diat_idch()
#' diat_ilm()
#' diat_lobo()
#' diat_sla()
#' diat_spear()
#' diat_pbidw()
#' diat_disp()
#' diat_idap()
#' diat_cemfgs_rb()
#' diat_checkName()
#' diat_getDiatBarcode()
#' diat_taxaList()
#'
#' @docType package
#' @name diaThor
#' @encoding UTF-8
#'

NULL
