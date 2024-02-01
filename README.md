# Analysis for Wepprich *et al.* (2024) manuscript on Ohio butterfly voltinism

This repository shares data and code to reproduce our analysis of changes in butterfly voltinism based on 1996-2022 statewide monitoring organized by the [Ohio Lepidopterists](https://www.ohiolepidopterists.org/).

**In the code directory**, follow the numbered R scripts to reproduce the analysis.
* 00_degree_days.R: Download daily weather data for each monitoring site
* 01_data_prep.R: Preparing monitoring data for analysis
* 02_fit_plot_gam.R: Fit seasonal phenological models using generalized additive models (GAM)
  * Fig. 2A and 2B
* 03_estimate_abundance_by_generation.R: Uses GAM to impute missing surveys, fits mixture models, and estimates abundance/phenology for each generation
  * Fig. 2C
* 04_prep_model_data.R: Creating variables for last generation size and weather covariates
  * Fig. 2D
  * Weather trends reported in Results
* 05_species_models.R: Main analysis of 30 species models (last generation size, overwinter growth rate, population trend)
  * Example figures for each species not included in manuscript
* 06_community_models.R: Main analysis of community models (last generation size, overwinter growth rate)
  * Tables S1 and S2
  * Fig. 3
  * Fig. S2 and S3
* 07_voltinism_traits.R: Analysis of species' traits correlated with population trends, phylogenetic analysis
  * Table S3
  * Fig. 4
  * Fig. S4 and S5
* 08_phenology_plot.R: Figure showing all species phenological and voltinism variation
  * Fig. S1

Unnumbered R scripts are for various utility functions 
* utils.R

**In the data directory**, we include raw data and intermediate products described with a metadata file. 
Model results from generalized additive models were too large for inclusion.


**In the figures directory**, we include the output from the analysis that we include in the manuscript.
Example figures from 05_species_models.R included here, but not in manuscript.
