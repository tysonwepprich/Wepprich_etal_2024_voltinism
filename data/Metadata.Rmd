# Metadata

**Raw data**
* data.trim.2023.csv: Every monitoring survey with site, date, species, count, duration of survey
* OHsites2023update.txt: Monitoring site geographic coordinates
* speciestraits.csv: Species traits and columns for inclusion in this analysis.
* species_names.csv: Used to match Latin names to common names
* OHbflyID.csv: Species list used in phylogenetic analysis
* RAxML_bipartitions.OHbflyBOLD_JRA_cds_aligned.gBlocks.consensus.famsConstrained.nex: Tree used in Wepprich et al. (2019), see paper's [supplement](https://doi.org/10.1371/journal.pone.0216270) for details.

**Intermediate products**
* daily_weather.rds: daily Daymet weather output from 00_degree_days.R, used in several other scripts.
* modparams.5.rds: output from 02_fit_plot_gam.R, such as fit and error terms for generalized additive models.
* mixmodcomparison.rds: output from 03_estimate_abundance_by_generation.R comparing mixture models
* genpops.rds: output from 03_estimate_abundance_by_generation.R with estimates for each species x site x year x generation
* covs.rds: weather covariates estimated in 04_prep_model_data.R for each species generation peak date
* modeling_data.rds: output of 04_prep_model_data.R, data with all covariates used in later scripts
* species_models.rds: output of 05_species_models.R containing parameters from the various species models
* pgls_bootstrap.rds: results of running PGLS 1000x with bootstrapped data in 07_voltinism_traits.R
* pgls_boot.csv: parameters from PGLS with bootstrap data
* pgls_noboot.csv: parameters from PGLS without data uncertainty accounted for.                            
                                                              

