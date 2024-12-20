
<!-- README.md is generated from README.Rmd. Please edit that file -->

# diaThor

*María Mercedes Nicolosi Gelis, Belén Sathicq, Joaquín Cochero*
<!-- badges: start --> <!-- badges: end -->

DiaThor calculates several diatom-based indices commonly used for water
quality assesment, directly from your species’ data

## Latest version: 0.1.4

We are bypassing the connection of diaThor to the DiatBarCode project
due to some URL changes in the latter. Until that is settled, this
version uses the internal DiatBarCode project’s database and will not
update to newer versions.

We also fixed some minor issues that made the diat_ips() function crash.

*Thanks to Julia Siegmund & Brent Bellinger (WPD, Austin, Texas) for
bringing this issues to us and helping us test the solution!*

## Description

The package calculates multiple biotic indices using diatoms from
environmental samples. Diatom species are recognized by their species’
name using a heuristic search, and their ecological data is retrieved
from multiple sources.

Morphological information about the species is retrieved from the
‘Diat.Barcode’ project:

- Rimet F., Gusev E., Kahlert M., Kelly M., Kulikovskiy M., Maltsev Y.,
  Mann D., Pfannkuchen M., Trobajo R., Vasselon V., Zimmermann J.,
  Bouchez A., 2019. Diat.barcode, an open-access curated barcode library
  for diatoms. Scientific Reports.
  <https://www.nature.com/articles/s41598-019-51500-6>

Size class classification is obtained from:

- Rimet F. & Bouchez A., 2012. Life-forms, cell-sizes and ecological
  guilds of diatoms in European rivers. Knowledge and management of
  aquatic ecosystems, 406: 1-14.
  <https://www.kmae-journal.org/articles/kmae/abs/2012/03/kmae120025/kmae120025.html>

Guild classification is obtained from:

- Rimet F. & Bouchez A., 2012. Life-forms, cell-sizes and ecological
  guilds of diatoms in European rivers. Knowledge and management of
  aquatic ecosystems, 406: 1-14.
  <https://www.kmae-journal.org/articles/kmae/abs/2012/03/kmae120025/kmae120025.html>

Ecological preferences are obtained form:

- Van Dam, H., Mertens, A., & Sinkeldam, J. (1994). A coded checklist
  and ecological indicator values of freshwater diatoms from the
  Netherlands. Netherland Journal of Aquatic Ecology, 28(1), 117-133.

Diversity index (Shannons H’) is calculated using the vegan package,
following:

- Shannon, C. E., and Weaver, W. (1949). ‘The Mathematical Theory of
  Communication.’ (University of Illinios Press: Urbana, IL, USA.)

Species tolerance and their ecological information to calculate each
biotic index is retrieved from their original sources:

- **IPS**: Coste, M. (1982). Étude des méthodes biologiques
  d’appréciation quantitative de la qualité des eaux. Rapport Cemagref
  QE Lyon-AF Bassin Rhône Méditerranée Corse.

- **TDI**: Kelly, M. G., & Whitton, B. A. (1995). The trophic diatom
  index: a new index for monitoring eutrophication in rivers. Journal of
  Applied Phycology, 7(4), 433-444.

- **IDP**: Gómez, N., & Licursi, M. (2001). The Pampean Diatom Index
  (IDP) for assessment of rivers and streams in Argentina. Aquatic
  Ecology, 35(2), 173-181.

- **DES**: Descy, J. P. 1979. A new approach to water quality estimation
  using diatom. Beih. Nov Hedw. 64:305-323

- **EPID**: Dell’Uomo, A. (1996). Assessment of water quality of an
  Apennine river as a pilot study for diatom-based monitoring of Italian
  watercourses. Use of algae for monitoring rivers, 65-72.

- **IDAP**: Prygiel, J., & Coste, M. (1993). The assessment of water
  quality in the Artois-Picardie water basin (France) by the use of
  diatom indices. Hydrobiologia, 269(1), 343-349.

- **ID-CH**: Hürlimann J., Niederhauser P. 2007: Méthodes d’analyse et
  d’appréciation des cours d’eau. Diatomées Niveau R (région). État de
  l’environnement n° 0740. Office fédéral de l’environnement, Berne. 132
  p

- **ILM**: Leclercq, L., & Maquet, B. (1987). Deux nouveaux indices
  diatomique et de qualité chimique des eaux courantes. Comparaison avec
  différents indices existants. Cahier de Biology Marine, 28, 303-310.

- **LOBO**:

  - Lobo, E. A., Callegaro, V. L. M., & Bender, E. P. (2002). Utilização
    de algas diatomáceas epilíticas como indicadoras da qualidade da
    água em rios e arroios da Região Hidrográfica do Guaíba, RS, Brasil.
    Edunisc.
  - Lobo, E. A., Bes, D., Tudesque, L., & Ector, L. (2004). Water
    quality assessment of the Pardinho River, RS, Brazil, using
    epilithic diatom assemblages and faecal coliforms as biological
    indicators. Vie et Milieu, 54(2-3), 115-126.

- **SLA**: Sládeček, V. (1986). Diatoms as indicators of organic
  pollution. Acta hydrochimica et hydrobiologica, 14(5), 555-566.

- **SPEAR(herbicides)**: Wood, R. J., Mitrovic, S. M., Lim, R. P.,
  Warne, M. S. J., Dunlop, J., & Kefford, B. J. (2019). Benthic diatoms
  as indicators of herbicide toxicity in rivers–A new SPEcies At Risk
  (SPEARherbicides) index. Ecological Indicators, 99, 203-213.

**Sample data included in the package is taken from**:

- Nicolosi Gelis, María Mercedes; Cochero, Joaquín; Donadelli, Jorge;
  Gómez, Nora. 2020. “Exploring the use of nuclear alterations, motility
  and ecological guilds in epipelic diatoms as biomonitoring tools for
  water quality improvement in urban impacted lowland streams”.
  Ecological Indicators, 110, 105951.

## Installation

You can install the released version of diathor from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("diathor")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("limnolab/DiaThor/")
```

## Example

To demonstrate the most common use of DiaThor, the package includes
sample data with the abundance of 164 diatom species in 108 sampled
sites (Nicolosi Gelis et al., 2020).

Install the package and load it into the R environment

``` r
> install.packages("diathor")
> library(diathor)
```

Load the internally included sample data

``` r
> data("diat_sampleData")
```

Run diaThorAll to get all the outputs from the sample data with the
default settings, and store the results into the “results” object, to
also retain the output within R

``` r
> results <- diaThorAll(diat_sampleData) #If the sample data was used
```

Note: The package will request the Input file an Output folder through a
dialog box

``` r
[1] "Select Input file"
[1] "Select Results folder"
```

After the Results folder is selected, all the calculations conducted
will be shown in the console

Optionally, run each individual function with the results of the
diat_loadData() function, for instance:

``` r
> loadedData <- diat_loadData() # load data with the diat_loadData() function
> results <- diat_ips(loadedData ) # use the diat_ips() function to calculate the IPS index with the loaded data
```

## CRAN vs. GitHub package

The package available in the CRAN repository does not automatically
update the internal diatom database from the ‘Diat.Barcode’ project. It
uses the internal version of such database, **which is currently v.10.1
published on 25-06-2021**

The GitHub version of the package will update to the most recent
database, if it is different to the current internal version.
