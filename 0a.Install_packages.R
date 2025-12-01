# install packages from Bioconductor and devtools
# to install phyloseq:

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

install.packages("devtools")
devtools::install_github("ChiLiubio/microeco")

library(devtools)
install_github("Russel88/MicEco")

library(devtools)
install_github("EWTekwa/Richness")

install.packages("devtools")
devtools::install_github("adw96/breakaway")
remotes::install_github("adw96/DivNet")

if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install("metagenomeSeq")

library(devtools)
devtools::install_github("vmikk/metagMisc")

library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

install.packages("remotes")
remotes::install_github("ctanes/adonisplus")