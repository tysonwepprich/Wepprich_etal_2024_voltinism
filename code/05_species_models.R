## Header ---------------------------
## Script name: 05_species_models.R
## Purpose of script: 
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## 
## ---


# 7.	Last generation relative size by photoperiod, latitude, temperature, precipitation

library(lme4)
library(merTools)
library(lmerTest) # lmer and step functions masked
library(ggeffects)
library(broom)
library(corrplot)
library(ggrepel)

source('01_data_prep.R')

moddat <- readRDS("modeling_data.rds")

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
         zyear = YearNum - 2009) %>% 
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


# Species models ----
outlist <- list()
sp <- sort(unique(moddat2$CommonName))
for (i in sp){
  print(i)
  tmp <- moddat2 %>% filter(CommonName == i) %>% 
    droplevels()
  
  mod <- lmer(lastloglambda ~ zdensann + zyear + (zordsite + zordann)^2 + 
                (1|SiteID),
              data = tmp, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  
  step_fm <- step(mod, alpha.fixed = .1)
  mod <- get_model(step_fm)
  print(summary(mod))
  
  # example figure
  print(try(ggpredict(mod, "zyear") %>% plot() + ggtitle(i)))
  yrplt <- ggpredict(mod, "zyear")
  plot(yrplt) + 
    scale_x_continuous(labels = c(1994,1999,2004,2009, 2014, 2019, 2024)) +
    labs(
      x = "Year",
      y = "Last generation size",
      title = "Voltinism change over time") +
    theme_bw(base_size = 24) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) 
  
  
  print(try(ggpredict(mod, c("zordsite", "zordann")) %>% plot() + ggtitle(i)))
  lgplt <- ggpredict(mod, c("zordann", "zordsite[-.21,.21]"))
  plot(lgplt) + 
    labs(
      x = "Penultimate generation peak date (annual variation)",
      y = "Last generation size",
      title = "Phenology changes voltinism") +
    theme_bw(base_size = 18) +
    scale_colour_brewer(palette = "Set1", labels = c("-1 SD (earlier)", "+1 SD (later)")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          legend.position = c(.2, .2)) + 
    guides(color=guide_legend(title="Mean site\npeak date"))
  
  
  
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
  
  tmp2 <- tmp %>% drop_na(zwinterann)
  tmp2 <- tmp2 %>% drop_na(ow_lam)
  
  # tmp <- tmp[-which(is.na(tmp$zwinterann)),]
  # tmp <- tmp[-which(is.na(tmp$ow_lam)),]
  
  ow2<- lmer(ow_lam ~ zdensann + zyear + (zlastsite + zlastann + zfrostann + zwinterann)^3 +
               # (1 + zdensann + zyear + | SiteID),
               # (1|Year) +
               (1|SiteID),
             data = tmp2, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  
  # ow2<- lmer(ow_lam ~ (zdensann + zyear + zlastsite + zlastann + zfrostann + zwinterann)^3 +
  #              # (1 + zdensann + zlastann + zlastsite | CommonName) +
  #              (1|SiteID),
  #            data = tmp, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  
  step_fm <- step(ow2, alpha.fixed = .1)
  ow2 <- get_model(step_fm)
  print(summary(ow2))
  print(try(ggpredict(ow2, c("zlastann", "zlastsite", "zwinterann")) %>% plot() + ggtitle(i)))
  # print(try(ggpredict(ow2, c("zlastann", "zwinterann", "zfrostann")) %>% plot() + ggtitle(i)))
  # print(try(ggpredict(ow2, c("zlastsite", "zfrostann")) %>% plot() + ggtitle(i)))
  
  owplt <- ggpredict(ow2, c("zlastann[-3:3]", "zlastsite[-.97,.97]", "zwinterann[-.94,.94]")) 
  # change attributes so facets labeled (3rd one only?)
  attr(owplt, "terms") <- c("Last generation size", "Site mean last generation size", "Winter temperature")
  
  pltdat1 <- plot(owplt) + 
    scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
    scale_fill_brewer(palette = "Set1", direction = -1) +
    # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
    labs(title = "Population growth rate ~ voltinism x winter", x = "Last generation size (annual variation)", y = "Overwinter growth rate") +
    theme_bw(base_size = 16) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    guides(color=guide_legend(title="Mean site\nlast generation\nsize"))
  
  
  
  
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
  
  # # Add van de pol within means centering models?
  
  fit0 <- lmer(lastloglambda ~ zdensann + zord + (1|SiteID), data = tmp)
  fit1 <- lmer(lastloglambda ~ zdensann + zordsite + zordann + (1|SiteID), data = tmp)
  fit2 <- lmer(lastloglambda ~ zdensann + zord + zordsite + (1|SiteID), data = tmp)
  
  vdp_param <- bind_rows(tidy(fit0) %>% mutate(model = "vdp_1"), 
                         tidy(fit1) %>% mutate(model = "vdp_2"),
                         tidy(fit2) %>% mutate(model = "vdp_3"))
  # print(vdp_param)
  
  
  tmp3 <- genpop %>%
    ungroup() %>%
    filter(CommonName == i) %>%
    filter(nyr >= 10, ObsSurvTotal >= 10, PosObsWeeks >= 3, trap.N.allgens >= 10, YearTotal >= 3) %>%
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>%
    mutate(zyear = YearNum - 2009,
           zmeanLL = MeanListLength,
           logtime = log(MeanDuration)) %>%
    filter(gen == 1) %>%
    mutate(Index = round(trap.N)) %>%
    droplevels()
  
  tmp3 <- left_join(tmp3, distinct(tmp[,c("SiteID", "zlastsite")]))
  
  modgen <- glmer(Index ~ zyear 
                  + zmeanLL +
                    (1 | SiteYear) +
                    (1 | SiteID) +
                    (1 | Year),
                  offset = logtime,
                  family = poisson(link = "log"),
                  data = tmp3,
                  control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
  
  print(summary(modgen))
  
  try(modgen2 <- glmer(Index ~ zyear * zlastsite #*region 
                       + zmeanLL +
                         (1 | SiteYear) +
                         (1 | SiteID) +
                         (1 | Year),
                       offset = logtime,
                       family = poisson(link = "log"),
                       data = tmp3,
                       control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7))))
  
  print(summary(modgen2))
  
  pop_param <- tidy(modgen) %>% mutate(model = "trend")
  popsite_param <- tidy(modgen2) %>% mutate(model = "sitetrend")
  
  # example trend plot
  tr2 = poptrend::ptrend(Index ~ trend(zyear, tempRE = TRUE, type = "loglinear") +
                           zmeanLL +
                           offset(logtime) +
                           s(SiteIDfact, bs = "re"),
                         family = nb(theta = NULL, link = "log"),
                         data = tmp3 %>% droplevels())
  
  plot(tr2)
  
  # analysis of predicted consequences of lost generations at sites with 10 or more yrs of observations
  species_lg_sd <- sd(tmp$zlastann, na.rm = TRUE)
  
  sitelist <- list()
  rm(aa)
  tmp2 <- tmp2 %>%
    group_by(SiteID) %>%
    mutate(nlams = length(ow_lam))
  sitereps <- unique(tmp2$SiteID[which(tmp2$nlams>=5)])
  
  
  for (site in sitereps){
    aa <- tmp2 %>% filter(SiteID == site) %>%
      arrange(Year) %>%
      dplyr::select(Year, SiteID, lastloglambda, zdensann, zyear, zlastsite, zlastann, zfrostann, zwinterann, ow_lam, nlams)
    
    # simulate owlam with different zlastann with model uncertainty
    # for each simulation compare observed, optimal, sitemean, and random zlastann, then geometric mean across years
    # output for each sim x site x scenario gm_mean, then compare differences (with uncertainty from sim) in scenarios for each site
    # for statewide, would it be a simple mean across sites?
    
    # PROBLEM: if zlastann not in the ow2 model, then it still selects "optimal" values that are just the 1st in the sequence.
    # FIX: don't do the comparison for species without zlastann in the model
    
    ## A function for simulating at new x-values
    simulateX <- function(object, nsim = 1, seed = NULL, X, ...) {
      object$fitted.values <- predict(object, X)
      simulate(object = object, nsim = nsim, seed = seed, ...)
    }
    
    
    # observed
    if(class(ow2) == "lm"){
      obssims <- simulateX(ow2,  nsim = 101, X = aa)
      # obspred <- predictInterval(ow2, newdata = aa, n.sims = 101, stat = "median", returnSims = TRUE)
      # obssims <- attr(obspred, "sim.results")
      
      # # optimal
      # # use wiggle to search for best zlastann
      # lg_range <- seq(from = min(aa$zlastann), to = max(aa$zlastann), length.out = 25)
      # range_pred <- wiggle(aa, varlist = "zlastann",
      #                      valueslist = list(lg_range))
      # range_pred$yhat <- predict(ow2, newdata = range_pred)
      # optannvar <- range_pred %>%
      #   group_by(Year) %>%
      #   summarise(best_lg = zlastann[which.max(yhat)])
      # 
      # # find max
      # # then use predictInterval
      # optdat <- aa %>%
      #   mutate(zlastann = optannvar$best_lg)
      # 
      # optsims <- simulateX(ow2,  nsim = 101, X = optdat)
      # 
      # mean_annual_diff = mean(aa$zlastann - optdat$zlastann) # negative means optimal LG is larger than observed
      # yrs_toosmall <- length(which((aa$zlastann - optdat$zlastann)<0))
      # yrs_toolarge <- length(which((aa$zlastann - optdat$zlastann)>0))
      
      
      # random
      randsims <- array(data = NA, dim = c(nrow(aa), 101))
      for (sim in 1:101){
        randpred <- aa %>%
          mutate(zlastann = sample(zlastann))
        randsim <- simulateX(ow2,  nsim = 1, X = randpred)
        randsims[,sim] <- randsim[,1]
      }
      
      # larger LG
      largedat <- aa %>%
        mutate(zlastann = zlastann + species_lg_sd)
      largesims <- simulateX(ow2,  nsim = 101, X = largedat)
      
      
      # smaller LG overall
      smalldat <- aa %>%
        mutate(zlastann = zlastann - species_lg_sd)
      smallsims <- simulateX(ow2,  nsim = 101, X = smalldat)
      
      
      # zlastann = 0
      meanpred <- aa %>%
        mutate(zlastann = 0)
      meansims <- simulateX(ow2,  nsim = 101, X = meanpred)
      
    }else{
      
      obspred <- predictInterval(ow2, newdata = aa, n.sims = 101, stat = "median", returnSims = TRUE)
      obssims <- attr(obspred, "sim.results")
      
      # # optimal
      # # use wiggle to search for best zlastann
      # lg_range <- seq(from = min(aa$zlastann), to = max(aa$zlastann), length.out = 25)
      # range_pred <- wiggle(aa, varlist = "zlastann",
      #                      valueslist = list(lg_range))
      # range_pred$yhat <- predict(ow2, newdata = range_pred)
      # optannvar <- range_pred %>%
      #   group_by(Year) %>%
      #   summarise(best_lg = zlastann[which.max(yhat)])
      # # find max
      # # then use predictInterval
      # optdat <- aa %>%
      #   mutate(zlastann = optannvar$best_lg)
      # optpred <- predictInterval(ow2, newdata = optdat, n.sims = 101, stat = "median", returnSims = TRUE)
      # optsims <- attr(optpred, "sim.results")
      # 
      # mean_annual_diff = mean(aa$zlastann - optdat$zlastann) # negative means optimal LG is larger than observed
      # yrs_toosmall <- length(which((aa$zlastann - optdat$zlastann)<0))
      # yrs_toolarge <- length(which((aa$zlastann - optdat$zlastann)>0))
      
      # larger LG
      largedat <- aa %>%
        mutate(zlastann = zlastann + species_lg_sd)
      largepred <- predictInterval(ow2, newdata = largedat, n.sims = 101, stat = "median", returnSims = TRUE)
      largesims <- attr(largepred, "sim.results")
      
      # smaller LG overall
      smalldat <- aa %>%
        mutate(zlastann = zlastann - species_lg_sd)
      smallpred <- predictInterval(ow2, newdata = smalldat, n.sims = 101, stat = "median", returnSims = TRUE)
      smallsims <- attr(smallpred, "sim.results")
      
      # # random
      # randsims <- array(data = NA, dim = c(nrow(aa), 101))
      # for (sim in 1:101){
      #   randdat <- aa %>%
      #     mutate(zlastann = sample(zlastann))
      #   randpred <- predictInterval(ow2, newdata = randdat, n.sims = 3, stat = "median", returnSims = TRUE)
      #   randsim <- attr(randpred, "sim.results")
      #   randsims[,sim] <- randsim[,1]
      # }
      #
      # randpred <- aa %>%
      #   mutate(zlastann = sample(zlastann))
      # randpred <- predictInterval(ow2, newdata = randpred, n.sims = 101, stat = "median", returnSims = TRUE)
      # randsims <- attr(randpred, "sim.results")
      #
      
      # zlastann = 0
      meandat <- aa %>%
        mutate(zlastann = 0)
      meanpred <- predictInterval(ow2, newdata = meandat, n.sims = 101, stat = "median", returnSims = TRUE)
      meansims <- attr(meanpred, "sim.results")
      
    } # close lm/lmer ifelse
    
    
    
    simdat <- rbind(rowMeans(log(apply(obssims, 2, FUN = function(x) Gmean(exp(x), conf.level = .5)))),
                    # rowMeans(log(apply(optsims, 2, FUN = function(x) Gmean(exp(x), conf.level = .5)))),
                    # rowMeans(log(apply(randsims, 2, FUN = function(x) Gmean(exp(x), conf.level = .5)))),
                    rowMeans(log(apply(largesims, 2, FUN = function(x) Gmean(exp(x), conf.level = .5)))),
                    rowMeans(log(apply(smallsims, 2, FUN = function(x) Gmean(exp(x), conf.level = .5)))),
                    rowMeans(log(apply(largesims - smallsims, 2, FUN = function(x) Gmean(exp(x), conf.level = .5)))),
                    rowMeans(log(apply(meansims, 2, FUN = function(x) Gmean(exp(x), conf.level = .5)))))
    simdat <- data.frame(simdat)
    simdat$scen <- c("obs", "large", "small", "diff", "mean")
    
    out <- data.frame(CommonName = i, SiteID = site, zlastsite = aa$zlastsite[1], nlams = aa$nlams[1],
                      simdat)
    
    sitelist[[(length(sitelist)+1)]] <- out
    
  }
  
  if(length(sitelist) > 0){
    site_lam <- bind_rows(sitelist)
    
    # note, could just use the "diff" scenario and it's uncertainty rather than below, but
    # I don't know how to do an lm with uncertainty on the response measurement which would be ideal to average across sites
    sitelg <- site_lam %>%
      group_by(SiteID, zlastsite, nlams) %>%
      summarise(largelam = mean[which(scen == "large")] - mean[which(scen == "obs")],
                smalllam = mean[which(scen == "small")] - mean[which(scen == "obs")])
    
    ggplot(sitelg, aes(x = zlastsite, y = largelam - smalllam))  + geom_point() + geom_smooth(method = "lm")
    
    sitelgmod1 <- lm(largelam - smalllam ~ 1, data = sitelg)
    site_param1 <- tidy(sitelgmod1) %>% mutate(model = "lg_sim1")
    
    sitelgmod <- lm(largelam - smalllam ~ zlastsite, data = sitelg)
    site_param <- tidy(sitelgmod) %>% mutate(model = "lg_sim2")
  }else{
    site_param <- data.frame(term = NA, estimate = NA, std.error = NA, statistic = NA, p.value = NA, model = "lg_sim")
  }
  outdf <- bind_rows(lastgen_param, vdp_param, ow_param, pop_param, popsite_param, site_param, site_param1) %>% mutate(CommonName = i)
  outlist[[(length(outlist)+1)]] <- outdf
  
}
allpars <- bind_rows(outlist)
saveRDS(allpars, "species_models.rds")

