# soil-biodiversity-food-forests
This repository contains all the code necessary to reproduce the analyses of soil physicochemical properties and soil biodiversity in temperate food forests as compared to forests, grasslands and croplands presented in the article "Planting food forests can increase soil biodiversity in agricultural landscapes of Northwest Europe".
The data used for these analyses are contained in the repository as well, but also archived on Zenodo at (https://10.5281/zenodo.17775912).

For the preparation of the sequencing data, please refer to the folder "Phyloseq_files", where you can find both the original data files (*.rds) and the code to process these into the files used for the analyses.

The folder "Data" contains the 
- information on the soil physicochemical properties (Chemphys.csv)
- raw data of the biomass (bacteria, non-AM fungi and AM-fungi: PLFA_NLFA), total group counts (nematodes, mites, springtails), and counts per species (macrofaunal groups)
- information on the samples and study sites (codes_info.csv, More_site_info.csv, pitfall_sampled.csv)
- the preprocessed *.rds files that contain estimates of Hill numbers per sample and Bray-Curtis dissimilarity matrices for the sequencing data
- information on the taxonomy of the macrofaunal groups

The main folder contains 
- the scripts with the code to perform the analyses for each of the four main hypotheses (H1, H2, H3 and H4) as well as for the analysis of the physicochemical soil properties; this also involves the file 0c.bootstrap_approach.R, which contains the code used for performing the bootstrapping, here just for one example.
- the files *.rds files produced and used over the course of the analyses (e.g. the bootstrap values in bblist_dens.rds)

The folder "Output_tables" contains tables with information on the fitted models and unstandardised effect sizes (ratios of means), which are produced using the scripts in the main folder.