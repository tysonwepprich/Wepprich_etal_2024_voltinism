# Analysis for Wepprich *et al.* (2024) manuscript on Ohio butterfly voltinism

This repository shares data and code to reproduce our analysis of changes in butterfly voltinism based on 1996-2022 statewide monitoring organized by the Ohio Lepidopterists.

**In the code directory**, follow the numbered R scripts to reproduce the analysis. 
* 00_degree_days.R 
* 01_monitoring_data_prep.R
* 02_fit_gam.R
* 03_plot_gam.R
* 04_fit_mixture_model.R
* 05_species_models.R
* 06_community_models.R
* 07_phylogenetic_models.R

Unnumbered R scripts are for utility functions and processing temperature rasters. 
* DDRP_cohorts_v1_funcs.R are functions for ...

**In the data directory**, we include raw data and intermediate products described with a metadata file. Model results from generalized additive models were too large for inclusion.

**In the figures directory**, we include the output from the analysis that we include in the manuscript.