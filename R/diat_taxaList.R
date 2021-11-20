#' Creates a single list with taxa names from all indices within DiaThor
#' @description
#' Creates a single list with taxa names from all indices within DiaThor
#' @keywords ecology diatom bioindicator biotic
#' @encoding UTF-8
#' @importFrom stringdist stringdist ain
#' @export diat_taxaList
#'


diat_taxaList <- function(){
  #gets internal indices
  ipsDB <- diathor::ips
  desDB <- diathor::des
  dispDB <- diathor::disp
  epidDB <- diathor::epid
  idapDB <- diathor::idap
  idchDB <- diathor::idch
  idpDB <- diathor::idp
  ilmDB <- diathor::ilm
  ipsDB <- diathor::ips
  loboDB <- diathor::lobo
  pbidwDB <- diathor::pbidw
  slaDB <- diathor::sla
  spearDB <- diathor::spear
  tdiDB <- diathor::tdi
  cemfgs_rbDB <- diathor::cemfgs_rb

  #create single vector list
  taxaList <- as.data.frame(c(ipsDB$fullspecies,
                              desDB$fullspecies,
                              dispDB$fullspecies,
                              epidDB$fullspecies,
                              idapDB$fullspecies,
                              idchDB$fullspecies,
                              idpDB$fullspecies,
                              ilmDB$fullspecies,
                              ipsDB$fullspecies,
                              loboDB$fullspecies,
                              pbidwDB$fullspecies,
                              slaDB$fullspecies,
                              spearDB$fullspecies,
                              tdiDB$fullspecies,
                              cemfgs_rbDB$fullspecies

  ))
  colnames(taxaList) <- "species"

  #trim blanks
  taxaList$species <- trimws(taxaList$species)
  #remove duplicates
  inner_taxaList <- as.data.frame(taxaList[!duplicated(taxaList[,"species"]),])
  colnames(inner_taxaList) <- "species"
  #order
  inner_taxaList <- as.data.frame(sort(inner_taxaList$species))
  colnames(inner_taxaList) <- "species"
  #return
  return(inner_taxaList)
}
