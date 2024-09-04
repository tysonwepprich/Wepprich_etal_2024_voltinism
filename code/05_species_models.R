## Header ---------------------------
## Script name: 05_species_models.R
## Purpose of script: Main analysis of 30 species models
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## Loop through models for each species and save results.
## Model of last generation size and its change in response to penultimate generation phenology
## Model of overwinter population growth rate in response to last generation size x winter onset/temperature
## Model of species statewide trends
## Simulation of population consequences of larger generations (integrating results from overwinter model)
## Optional plots of species results (not included in manuscript)
## Results saved in "data/species_models.rds", which is used for the voltinism trait vs long-term trend analysis in 07_voltinism_traits.R.
## Alternative set of results run with imputation weighting in models, checked differences at bottom of this script
## ---

library(lme4)
library(merTools)
library(lmerTest) # lmer and step functions masked
library(ggeffects)
library(broom)
library(broom.mixed)
library(DescTools)
library(performance)


source('code/01_data_prep.R')

moddat <- readRDS("data/modeling_data.rds")
genpop <- readRDS("data/genpops.rds")

# Within-group centering ----
# predictors as species group-level, site group-level, and annual variation
moddat2 <- moddat %>% 
  ungroup() %>% 
  mutate(SpSiteYr = paste(CommonName, SiteYear, sep = "_"),
         response = round(lastN),
         response1 = round(firstN),
         ow_lam = log((firstN + 1)/(trap.N + 1)),
         off = log(round(trap.N)),
         zord = scale_this(ord),
         zphoto = scale_this(photo),
         zfrost = scale_this(firstfrostdoy),
         zwinter = scale_this(allwintertemp),
         zyear = YearNum - 2009,
         # Covariates added for imputation to test if imputed surveys are unduly affecting results
         last_imputation_weights = (lastObs + 1)/(lastImp + lastObs + 1), # NA's if both Imp/Obs are zero
         first_imputation_weights = (firstObs + 1)/(firstImp + firstObs + 1)) %>% # NA's if both Imp/Obs are zero
  group_by(CommonName) %>% 
  mutate(zordspec = mean(zord),
         zphotospec = mean(zphoto),
         zdenspec = mean(off),
         zlastspec = mean(lastloglambda)) %>% 
  group_by(SiteID, CommonName) %>% 
  mutate(
    zordsite = mean(zord) - zordspec,
    zordann = zord - mean(zord),
    zphotosite = mean(zphoto) - zphotospec,
    zphotoann = zphoto - mean(zphoto),
    zdensite = mean(off) - zdenspec,
    zdensann = off - mean(off),
    zlastsite = mean(lastloglambda) - zlastspec,
    zlastann = lastloglambda - mean(lastloglambda),
    nlams = length(which(!is.na(ow_lam)))) %>% 
  group_by(SiteID) %>% 
  mutate(zfrostsite = mean(zfrost),
         zfrostann = zfrost - mean(zfrost),
         zwintersite = mean(zwinter, na.rm = TRUE),
         zwinterann = zwinter - mean(zwinter, na.rm = TRUE))

latin <- read.csv("data/species_names.csv") %>% 
  dplyr::select(CommonName, Genus, Species) %>%
  mutate(Latin = paste(Genus, Species, sep = " "))

moddat2 <- left_join(moddat2, latin)

# Imputation weights by species based on observed counts versus imputed counts
moddat2 %>% group_by(CommonName) %>% 
  summarise(last = sum(response, na.rm = TRUE), 
            last_obs = sum(response * last_imputation_weights), 
            last_perc_obs = last_obs / (last),
            first = sum(response1, na.rm = TRUE),
            first_obs = sum(response1 * first_imputation_weights, na.rm = TRUE),
            first_perc_obs = first_obs / first) %>% 
  summary()

# Fig S13: Imputation weights ----
hist(moddat2$first_imputation_weights, main = "First generation imputation", xlab = "Proportion of count observed vs. imputed")
hist(moddat2$last_imputation_weights, main = "Last generation imputation", xlab = "Proportion of count observed vs. imputed")

# Species models ----
outlist <- list()
sp <- sort(unique(moddat2$CommonName))
for (i in sp){
  print(i)
  tmp <- moddat2 %>% filter(CommonName == i) %>% 
    droplevels()
  
  latin <- tmp$Latin[1]
  
  # Last generation size model ----
  mod <- lmer(lastloglambda ~ zdensann + zyear + (zordsite + zordann)^2 + 
                (1|SiteID),
              data = tmp, 
              # weights = last_imputation_weights,
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  
  step_fm <- step(mod, alpha.fixed = .1, keep = c("zdensann", "zyear"))
  mod <- get_model(step_fm)
  # print(summary(mod))
  # check_model(mod)
  
  # # example figures
  # print(try(ggpredict(mod, "zyear") %>% plot() + ggtitle(i)))
  # yrplt <- ggpredict(mod, "zyear")
  # yrplt <- plot(yrplt, add.data = TRUE) +
  #   scale_x_continuous(breaks = c(-9, -4, 1, 6, 11), labels = c(2000, 2005, 2010, 2015, 2020), limits = c(-14,14), expand = c(0,0)) +
  #   labs(
  #     x = "Year",
  #     y = "Last generation size",
  #     title = NULL) +
  #   theme_bw(base_size = 18) +
  #   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  # yrplt
  # saveRDS(yrplt, file = paste0("figures/species_examples/lgtrend_", latin, ".rds"))
  

  # print(try(ggpredict(mod, c("zordsite", "zordann")) %>% plot() + ggtitle(i)))
  # lgplt <- ggpredict(mod, c("zordann", "zordsite[-.22,.22]"))
  # lgplt <- plot(lgplt) +
  #   labs(
  #     x = "Penultimate generation peak date (annual variation SD)",
  #     y = "Last generation size",
  #     title = NULL) +
  #   theme_bw(base_size = 18) +    
  #   # scale_colour_manual(values = c("#E41A1C", "#377EB8"), labels = c("-1 SD (earlier)", "+1 SD (later)")) +
  #   scale_colour_brewer(palette = "Set1", direction = 1, labels = c("-1 SD (earlier)", "+1 SD (later)")) +
  #   scale_fill_brewer(palette = "Set1", direction = 1, labels = c("-1 SD (earlier)", "+1 SD (later)")) +
  #   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
  #         legend.position = c(.2, .2)) +
  #   guides(color=guide_legend(title="Mean site\npeak date"))
  # lgplt
  # saveRDS(lgplt, file = paste0("figures/species_examples/lgresp_", latin, ".rds"))
  
  
  
  if (class(mod) == "lmerModLmerTest"){
    r2 <- MuMIn::r.squaredGLMM(mod)
    lastgen_param <- tidy(mod) %>% 
      mutate(model = "lastgen",
             # overdisp = blmeco::dispersion_glmer(mod),
             marg_r2 = r2[1,1],
             cond_r2 = r2[1,2])
  }else{
    r2 <- glance(mod)$r.squared
    lastgen_param <- tidy(mod) %>% 
      mutate(model = "lastgen",
             # overdisp = NA,
             marg_r2 = r2,
             cond_r2 = NA)
  }
  
  
  
  # Overwinter population growth rate model ----
  tmp2 <- tmp %>% drop_na(zwinterann)
  tmp2 <- tmp2 %>% drop_na(ow_lam)
  
  ow2<- lmer(ow_lam ~ zdensann + zyear + (zlastsite + zlastann + zfrostann + zwinterann)^3 +
               # (1 + zdensann + zyear + | SiteID),
               # (1|Year) +
               (1|SiteID),
             # weights = last_imputation_weights + first_imputation_weights,
             data = tmp2, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  
  step_fm <- step(ow2, alpha.fixed = .1)
  ow2 <- get_model(step_fm)
  # print(summary(ow2))
  # check_model(ow2)
  
  # # Example figure
  # # print(try(ggpredict(ow2, c("zlastann", "zlastsite", "zwinterann")) %>% plot() + ggtitle(i)))
  # 
  # owplt <- ggpredict(ow2, c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-1,1]"))
  # # change attributes so facets labeled (3rd one only?)
  # attr(owplt, "terms") <- c("Last generation size", "Site mean last generation size", "Winter temperature")
  # 
  # owplt <- plot(owplt) +
  #   scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  #   scale_fill_brewer(palette = "Set1", direction = -1) +
  #   # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  #   labs(title = NULL, x = "Last generation size (annual variation SD)", y = "Overwinter growth rate") +
  #   theme_bw(base_size = 18) +
  #   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  #   theme(legend.position = c(.125, .15)) +
  #   guides(color=guide_legend(title="Mean site last\ngeneration size"))
  # owplt
  # saveRDS(owplt, file = paste0("figures/species_examples/ow_", latin, ".rds"))
  
  
  if (class(mod) == "lmerModLmerTest"){
    r2 <- MuMIn::r.squaredGLMM(ow2)
    ow_param <- tidy(ow2) %>% 
      mutate(model = "ow",
             # overdisp = blmeco::dispersion_glmer(mod),
             marg_r2 = r2[1,1],
             cond_r2 = r2[1,2])
  }else{
    r2 <- glance(ow2)$r.squared
    ow_param <- tidy(ow2) %>% 
      mutate(model = "ow",
             # overdisp = NA,
             marg_r2 = r2,
             cond_r2 = NA)
  }
  
  # Local adaptation parameters ----
  # van de pol within-means centering models
  # response of last generation size to spatial means and annual variation
  
  fit0 <- lmer(lastloglambda ~ zdensann + zord + (1|SiteID), data = tmp)
  # equation 2 from van de pol & wright (2009)
  # coef of zordsite == between site sensitivity, zordann == within site/between year sensitivity
  fit1 <- lmer(lastloglambda ~ zdensann + zordsite + zordann + (1|SiteID), data = tmp)
  # equation 3 from van de pol & wright (2009)
  # coefficient for zordsite == Between - within group difference
  # (equivalent to space - time sensitivity in Ramirez-Pareda et al (2024))
  fit2 <- lmer(lastloglambda ~ zdensann + zord + zordsite + (1|SiteID), data = tmp)

  # # with imputation weights
  # fit0 <- lmer(lastloglambda ~ zdensann + zord + (1|SiteID), data = tmp, weights = last_imputation_weights)
  # fit1 <- lmer(lastloglambda ~ zdensann + zordsite + zordann + (1|SiteID), data = tmp, weights = last_imputation_weights)
  # fit2 <- lmer(lastloglambda ~ zdensann + zord + zordsite + (1|SiteID), data = tmp, weights = last_imputation_weights)

  vdp_param <- bind_rows(tidy(fit0) %>% mutate(model = "vdp_1"), 
                         tidy(fit1) %>% mutate(model = "vdp_2"),
                         tidy(fit2) %>% mutate(model = "vdp_3"))
  # print(vdp_param)
  
  
  # Species' statewide population trend ----
  tmp3 <- genpop %>%
    ungroup() %>%
    filter(CommonName == i) %>%
    filter(nyr >= 10, ObsSurvTotal >= 10, PosObsWeeks >= 3, trap.N.allgens >= 10, YearTotal >= 3) %>%
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>%
    mutate(zyear = YearNum - 2009,
           zmeanLL = MeanListLength,
           logtime = log(MeanDuration)) %>%
    filter(gen == 1) %>%
    mutate(Index = round(trap.N),
           first_imputation_weights = (total.obs.gen + 1)/(total.imp.gen + total.obs.gen + 1)) %>%
    droplevels()
  
  modgen <- glmer(Index ~ zyear 
                  + zmeanLL +
                    (1 | SiteYear) +
                    (1 | SiteID) +
                    (1 | Year),
                  offset = logtime,
                  family = poisson(link = "log"),
                  data = tmp3,
                  # weights = first_imputation_weights,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
  
  # print(summary(modgen))
  pop_param <- tidy(modgen) %>% mutate(model = "trend")
  
  
  # # example trend plot with annual variation (random effects)
  # newdat <- data.frame(zyear = sort(unique(tmp3$zyear)), Year = as.character(c(1996:2022)),
  #                      zmeanLL = mean(tmp3$zmeanLL), logtime = mean(tmp3$logtime))
  # newdat$annvar <- predict(modgen, newdat, re.form = ~ (1 | Year), type = "response")
  # 
  
  # print(try(ggpredict(modgen, c("zyear")) %>% plot() + ggtitle(i)))
  
  # trendplt <- ggpredict(modgen, c("zyear"))
  
  # trendplt <- plot(trendplt) +
  #   scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  #   scale_fill_brewer(palette = "Set1", direction = -1) +
  #   geom_point(data = newdat, aes(x = zyear, y = annvar)) +
  #   scale_x_continuous(breaks = c(-9, -4, 1, 6, 11), labels = c(2000, 2005, 2010, 2015, 2020), limits = c(-14,14), expand = c(0,0)) +
  #   scale_y_continuous(limits=c(0,NA), expand = expansion(mult = c(0, .1))) +
  #   labs(title = NULL, x = "Year", y = "Population Index (1st generation)") +
  #   theme_bw(base_size = 18) +
  #   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  # trendplt
  # 
  # saveRDS(trendplt, file = paste0("figures/species_examples/trend_", latin, ".rds"))
  
  # Simulation of consequences of larger last generation ----
  species_lg_sd <- sd(tmp$zlastann, na.rm = TRUE)
  
    # simulate owlam with different zlastann with model uncertainty
    # for each simulation compare observed, sitemean, and random zlastann, then geometric mean across years
    # output for each sim x site x scenario gm_mean, then compare differences (with uncertainty from sim) in scenarios for each site
    
    ## A function for simulating at new x-values
    simulateX <- function(object, nsim = 1, seed = NULL, X, ...) {
      object$fitted.values <- predict(object, X)
      simulate(object = object, nsim = nsim, seed = seed, ...)
    }
    
    # observed  # try tmp2 (not aa)
    if(class(ow2) == "lm"){
      obssims <- simulateX(ow2,  nsim = 101, X = tmp2)
      
      # larger LG
      largedat <- tmp2 %>%
        mutate(zlastann = zlastann + species_lg_sd)
      largesims <- simulateX(ow2,  nsim = 101, X = largedat)
      
      
      # smaller LG overall
      smalldat <- tmp2 %>%
        mutate(zlastann = zlastann - species_lg_sd)
      smallsims <- simulateX(ow2,  nsim = 101, X = smalldat)
      
      
      # zlastann = 0
      meanpred <- tmp2 %>%
        mutate(zlastann = 0)
      meansims <- simulateX(ow2,  nsim = 101, X = meanpred)
      
    }else{
      
      obspred <- predictInterval(ow2, newdata = tmp2, n.sims = 101, stat = "median", returnSims = TRUE)
      obssims <- attr(obspred, "sim.results")
      
      # larger LG
      largedat <- tmp2 %>%
        mutate(zlastann = zlastann + species_lg_sd)
      largepred <- predictInterval(ow2, newdata = largedat, n.sims = 101, stat = "median", returnSims = TRUE)
      largesims <- attr(largepred, "sim.results")
      
      # smaller LG overall
      smalldat <- tmp2 %>%
        mutate(zlastann = zlastann - species_lg_sd)
      smallpred <- predictInterval(ow2, newdata = smalldat, n.sims = 101, stat = "median", returnSims = TRUE)
      smallsims <- attr(smallpred, "sim.results")
      
      # zlastann = 0
      meandat <- tmp2 %>%
        mutate(zlastann = 0)
      meanpred <- predictInterval(ow2, newdata = meandat, n.sims = 101, stat = "median", returnSims = TRUE)
      meansims <- attr(meanpred, "sim.results")
      
    } # close lm/lmer ifelse
    
    
    
    simdat <- rbind(rowMeans(log(apply(obssims, 2, FUN = function(x) Gmean(exp(x), conf.level = .95)))),
                    rowMeans(log(apply(largesims, 2, FUN = function(x) Gmean(exp(x), conf.level = .95)))),
                    rowMeans(log(apply(smallsims, 2, FUN = function(x) Gmean(exp(x), conf.level = .95)))),
                    rowMeans(log(apply(largesims - smallsims, 2, FUN = function(x) Gmean(exp(x), conf.level = .95)))),
                    rowMeans(log(apply(meansims, 2, FUN = function(x) Gmean(exp(x), conf.level = .95)))))
    simdat <- data.frame(simdat)
    simdat$scen <- c("obs", "large", "small", "diff", "mean")
    
    out <- simdat %>% filter(scen == "diff")
    site_param <- data.frame(term = "diff", 
                             estimate = out$mean,
                             std.error = (out$upr.ci - out$mean) / 1.959964, 
                             statistic = NA, p.value = NA, model = "lg_sim1")
    
  outdf <- bind_rows(lastgen_param, vdp_param, ow_param, pop_param, site_param) %>% mutate(CommonName = i)
  outlist[[(length(outlist)+1)]] <- outdf
  
}
allpars <- bind_rows(outlist)
saveRDS(allpars, "data/species_models.rds")
# saveRDS(allpars, "data/species_models_check_impute.rds")


# Imputation robustness check ----
# Do parameter estimates change when imputed weekly surveys are weighted lower in the models?

orig_analysis <- readRDS("data/species_models.rds")
impute_analysis <- readRDS("data/species_models_check_impute.rds")

org <- orig_analysis[,c("term", "estimate", "model", "CommonName")] %>% 
  mutate(org_est = estimate) %>% 
  dplyr::select(-estimate) %>% 
  filter(term != "sd__(Intercept)", term != "sd__Observation")
imp <- impute_analysis %>% 
  mutate(imp_est = estimate) %>% 
  dplyr::select("term", "imp_est", "model", "CommonName") %>% 
  filter(term != "sd__(Intercept)", term != "sd__Observation")

test <- full_join(org, imp)

# 37 rows with NAs where imputed parameter estimates are not made (due to model selection)
test %>% 
  filter(model %in% c("lastgen", "ow")) %>% 
  group_by(term, model) %>% 
  summarise(r = cor(org_est, imp_est, use = "complete.obs"),
            missing_org = length(which(is.na(org_est))),
            missing_imp = length(which(is.na(imp_est)))) %>% 
  data.frame()

# check out individual model terms across species
test2 <- test %>% filter(term == "zyear", model == "ow")

ggplot(test2, aes(org_est, imp_est, label = CommonName)) +
  geom_point() +
  ggrepel::geom_text_repel() +
  geom_abline(slope = 1, intercept = 0)
