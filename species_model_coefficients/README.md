These three files contain results that did not fit in the Supporting Information.

* coefficients_last_gen_size.csv: Model results for all 30 species models of last generation size (voltinism shift). Output was produced from linear mixed effects models using the broom R package.
  Model terms:
  * (Intercept): model intercept
  * zdensann: density-dependence (scaled deviation from site's mean population size in penultimate generation)
  * zyear: temporal trend in last generation size
  * zordsite: mean site phenology of penultimate generation
  * zordann: annual anomaly in penultimate generation phenology (compared to site mean)
  * zordsite:zordann: interaction between the above two terms
  * sd_(Intercept): standard deviation of random intercept for each Site
  * sd_(Observation): standard deviation of residual variation
  
* coefficients_overwinter.csv: Model results for all 30 species models of overwinter population growth rates. Output was produced from linear mixed effects models using the broom R package.
  Model terms:
  * (Intercept): model intercept
  * zdensann: density-dependence (scaled deviation from site's mean population size in penultimate generation)
  * zyear: temporal trend in overwinter population growth rate
  * zlastsite: mean site last generation size
  * zlastann: annual anomaly in last generation size (compared to site mean)
  * zfrostann: annual anomaly in winter onset (first hard frost <-2C)
  * zwinterann: annual anomaly in mean minimum winter temperature (Nov. 1 - Mar. 31)
  * Plus 2 and 3-way interactions between the above terms
  * sd_(Intercept): standard deviation of random intercept for each Site
  * sd_(Observation): standard deviation of residual variation

* TraitsTable.csv: Voltinism traits that correspond to estimates and 95% confidence intervals in Fig. 3.
 * Latin: species name
 * Sample: data used for the species models of last generation size and overwinter population growth rates respectively.
 * Generations: statewide distribution of counts to different generations by percent.
 * meandoy: species' mean date of penultimate generation phenology
 * LG_trend: temporal trend in last generation size
 * LGspatial: sensitivity of last generation size to penultimate generation phenology consistent with local adaptation
 * LGannual: sensitivity of last generation size to penultimate generation phenology consistent with phenotypic plasticity
 * LGsim: geometric mean of simulated overwinter population growth rates if species had larger last generations across all observations
 * PopTrend: long-term statewide population trend estimated from the first generations only
