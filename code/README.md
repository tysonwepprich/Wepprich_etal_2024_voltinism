Follow the numbered R scripts to reproduce the analysis.

* 00_degree_days.R: Download daily weather data for each monitoring site
* 01_data_prep.R: Preparing monitoring data for analysis
  * Fig. S12
* 02_fit_plot_gam.R: Fit seasonal phenological models using generalized additive models (GAM)
  * Fig. 2A and 2B
  * All species phenology figures
* 03_estimate_abundance_by_generation.R: Uses GAM to impute missing surveys, fits mixture models, and estimates abundance/phenology for each generation
  * Fig. 2C
  * All species mixture model figures
* 04_prep_model_data.R: Creating variables for last generation size and weather covariates
  * Fig. 2D
  * Weather trends reported in Results
* 05_species_models.R: Main analysis of 30 species models (last generation size, overwinter growth rate, population trend)
  * Fig. S13
* 06_community_models.R: Main analysis of community models (last generation size, overwinter growth rate)
  * Tables S1 and S2 for model coefficients
  * Fig. 4
  * Fig. S3, S4, and S5
  * Tables S6 and S7 for imputation check
* 07_voltinism_traits.R: Analysis of species' traits correlated with population trends, phylogenetic analysis
  * Table S3
  * Fig. 3
  * Fig. 5A, B, C
  * Fig. S6, S7, S8, S9
  * Coefficients from species models
* 08_phenology_plot.R: Figure showing all species phenological and voltinism variation
  * Fig. S2
* 09_sitemap.R: Makes map of monitoring sites
  * Fig. S1

Unnumbered R scripts

* utils.R
  * Various utility functions 
* small_copper_example.R
  * Used in Supporting Information to describe mixture model assumptions
