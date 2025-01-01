## Header ---------------------------
## Script name: 03_estimate_abundance_by_generation.R
## Purpose of script: Uses GAM to impute missing surveys, fits mixture models, and estimates abundance/phenology for each generation
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## Running this script is not necessary to reproduce the main results, because
## output of this script is saved as "data/genpops.rds", which is used in subsequent scripts.
## You would need to run 02_fit_plot_gam.R to get required input for these mixture models
## ---

source('code/01_data_prep.R')

# Candidate species ----
# 42 candidate species for mixture models, not all will work (mv_analysis == 1)
# After species filter (mixture models work/sample size sufficient) 30 remain (final_mv_analysis == 1)
mvspecies <- read.csv("data/speciestraits.csv", header = TRUE) %>%
  # filter(mv_analysis == 1) %>%
  filter(final_mv_analysis == 1) %>%
  droplevels.data.frame()

modfiles <- list.files("gams", full.names = TRUE)
modfiles <- modfiles[grep(pattern = "final", x = modfiles, fixed = TRUE)]

splist <- list()
datlist <- list()
for (sp in 1:nrow(mvspecies)){
  species <- mvspecies$CommonName[sp]
  minvolt <- 1
  maxvolt <- mvspecies$Voltinism[sp]
  
  # Import GAM phenology model ---- 
  # and the count data used to fit it
  mf <- modfiles[grep(pattern = species, x = modfiles, fixed = TRUE)]
  tmp <- readRDS(mf)
  mod <- tmp$gammod
  dat <- tmp$datGAM %>% 
    filter(DOY > 80 & DOY < 315) # removing anchoring zeros used in GAM fitting
  
  # Impute missing surveys ----
  surv <- dat %>%
    dplyr::select(SiteID, Week, Year, SiteYear) %>%
    distinct()
  allsurv <- surv %>%
    complete(SiteID, Week, Year) %>%
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>%
    filter(Week <= 30) %>%
    filter(SiteYear %in% unique(surv$SiteYear))
  missing <- anti_join(allsurv, surv)

  impute_days <- expand.grid(Year = unique(missing$Year), Week = c(1:30)) %>%
    mutate(DOY = ifelse(as.numeric(as.character(Year)) %% 4 == 0, 92 + Week * 7 - 4, 91 + Week * 7 - 4)) # midweek to impute

  # Predict missing weeks.
  missing <- missing %>%
    left_join(impute_days) %>%
    mutate(SiteDate = as.Date(paste(DOY, Year, sep = "-"), format = "%j-%Y")) %>%
    left_join(gdd[, c("SiteID", "SiteDate", "lat", "lon", "accumdegday", "region")], by = c("SiteID", "SiteDate")) %>%
    mutate(RegYear = paste(region, Year, sep = "_")) %>%
    mutate(Total = predict(mod, newdata = ., type = "response"),
           datatype = "imputed")
  
  alldat <- dat %>%
    dplyr::select(SiteID, Week, Year, SiteYear, DOY, SiteDate, lat, lon, accumdegday, listlength, duration, region, RegYear, Total) %>% 
    mutate(datatype = "observed") %>%
    bind_rows(missing) %>%
    mutate(CommonName = species)
  
  # Mixture models with mclust package
  # compare 3 timescales: degree day, ordinal day, and degree day with log(counts)
  dd_gdd <- rep(alldat$accumdegday, round(alldat$Total))
  dd_doy <- rep(alldat$DOY, round(alldat$Total))
  dd_gddlog <- rep(alldat$accumdegday, round(log(alldat$Total + 1)))
  dd_doylog <- rep(alldat$DOY, round(log(alldat$Total + 1)))
  
  mm_gdd <- Mclust(dd_gdd, G = minvolt:maxvolt, modelNames = "E")
  mm_doy <- Mclust(dd_doy, G = minvolt:maxvolt, modelNames = "E")
  mm_gddlog <- Mclust(dd_gddlog, G = minvolt:maxvolt, modelNames = "E")
  mm_doylog <- Mclust(dd_doylog, G = minvolt:maxvolt, modelNames = "E")
  
  # test doy vs gdd
  # constraining to equal variance here, but could try unequal for more realism
  mmdf1 <- data.frame(species = species, gen = 1:length(mm_gdd$parameters$mean),
                     prop = mm_gdd$parameters$pro, mean = mm_gdd$parameters$mean,
                     time = "gdd", uncertainty = mean(mm_gdd$uncertainty))
  mmdf2 <- data.frame(species = species, gen = 1:length(mm_doy$parameters$mean),
                      prop = mm_doy$parameters$pro, mean = mm_doy$parameters$mean,
                      time = "doy", uncertainty = mean(mm_doy$uncertainty))
  mmdf3 <- data.frame(species = species, gen = 1:length(mm_gddlog$parameters$mean),
                      prop = mm_gddlog$parameters$pro, mean = mm_gddlog$parameters$mean,
                      time = "gdd_log", uncertainty = mean(mm_gddlog$uncertainty))
  mmdf4 <- data.frame(species = species, gen = 1:length(mm_doylog$parameters$mean),
                      prop = mm_doylog$parameters$pro, mean = mm_doylog$parameters$mean,
                      time = "doy_log", uncertainty = mean(mm_doylog$uncertainty))
  splist[[sp]] <- bind_rows(mmdf1, mmdf2, mmdf3, mmdf4)
  
  # check predictions
  mm <- list(mm_gdd, mm_doy, mm_gddlog, mm_doylog) 
  mm_results <- c(mmdf1$uncertainty[1], mmdf2$uncertainty[1], mmdf3$uncertainty[1], mmdf4$uncertainty[1])
  mm_results[which(mm_results == 0)] <- NA # only 1 generation means zero uncertainty
  mm <- mm[base::order(mm_results, na.last = TRUE)][[1]] # pick best model
  
 # Assigning generation ----
 # avoiding rounding so GAM predictions are multiplied by mixed model generation probabilities
 # each count split
  dat <- alldat %>% 
    group_by(SiteYear) %>% 
    mutate(YearTotal = sum(Total[which(datatype == "observed")]),
           ObsSurvTotal = length(Week[which(datatype == "observed")]),
           ImpSurvTotal = length(Week[which(datatype == "imputed")]),
           PosObsWeeks = length(Week[which(datatype == "observed" & Total > 0)]),
           MeanListLength = mean(listlength, na.rm = TRUE),
           MeanDuration = mean(duration, na.rm = TRUE))
  
  if(base::order(mm_results, na.last = TRUE)[1] %in% c(1,3)){
    mmpred <- predict(mm, newdata = dat$accumdegday)
  }else{
    mmpred <- predict(mm, newdata = dat$DOY)
  }  
  res <- as.data.frame(mmpred$z)
  
  dat2 <- cbind(dat, dat$Total * res)
  dat3 <- dat2 %>% 
    pivot_longer(cols = as.character(c(1:length(mm$parameters$mean))),
                 names_to = "gen",
                 values_to = "GenTotal")
  
  dat3$GenTotal <- round(dat3$GenTotal, digits = 3)
  
  
  # # Plot mixture model for all observations ----
  # mmplt <- ggplot(dat3 %>% filter(GenTotal > 0), aes(x = accumdegday, y = GenTotal, color = as.factor(gen), , shape = datatype)) +
  #   geom_point(alpha = .2) +
  #   # scale_y_continuous(limits = c(0, 190), expand = c(0,0), breaks = c(0, 50, 100, 150)) +
  #   # scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  #   scale_color_brewer(name = NULL, type = "qual", palette = "Dark2") +
  #   labs(color = "Generation") +
  #   labs(x = "Degree-days accumulated (5/30°C thresholds)") +
  #   labs(y = "Observed or imputed weekly count") +
  #   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  #   theme(legend.position = "none")
  # # mmplt
  
  colordf <- data.frame(gen = 1:length(mm$parameters$mean), mean = mm$parameters$mean, colorhex = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")[1:length(mm$parameters$mean)])
  
  # # Plot mixture model for by year for one region ----
  # # Example for supplement
  # mmplt <- ggplot(dat3 %>% filter(GenTotal >= 1, region == "NE", Year %in% c(2008:2017)), aes(x = accumdegday, y = GenTotal, color = as.factor(gen))) +
  #   geom_point(alpha = .5) +
  #   # scale_y_continuous(limits = c(0, 190), expand = c(0,0), breaks = c(0, 50, 100, 150)) +
  #   scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
  #   scale_color_brewer(name = NULL, type = "qual", palette = "Dark2") +
  #   geom_vline(data = colordf, aes(xintercept = mean), color = "black", linetype = "dashed") +
  #   labs(color = "Generation") +
  #   labs(x = "Degree-days accumulated (5/30°C thresholds)") +
  #   labs(y = "Observed or imputed weekly count") +
  #   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  #   theme(legend.position = "none") +
  #   facet_wrap(. ~ Year, ncol = 2, scales = "free_y")
  # # mmplt
  
  
 # Plot mixture model results by region and year ----  
  # name for figures
  latin <- read.csv("data/species_names.csv") %>% 
    dplyr::select(CommonName, Genus, Species) %>%
    filter(CommonName == tmp$params$CommonName) %>% 
    mutate(Latin = paste(Genus, Species, sep = " ")) %>% 
    pull(Latin) %>% 
    as.character()

  # Fig. 2C ----
  dat4 <- dat3 %>%
    ungroup() %>%
    filter(GenTotal >= 1, Year %in% c(2008:2017)) %>%
    mutate(yrgroup = ifelse(Year %in% c(2010, 2011, 2012, 2016, 2017), "Warm years", "Cool years"),
           sitegroup = ifelse(region %in% c("NE", "NW"), "Cool sites", "Warm sites"))

  mmplt <- ggplot(dat4, aes(x = accumdegday, y = GenTotal, color = as.factor(gen))) +
    geom_point(alpha = .3, size = .9) +
    scale_y_continuous(limits = c(0, 190), expand = c(0,0), breaks = c(0, 50, 100, 150)) +
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) +
    scale_color_brewer(name = NULL, type = "qual", palette = "Dark2") +
    facet_grid(sitegroup~yrgroup, scales = "free_y") +
    # facet_wrap(sitegroup~yrgroup, scales = "free_y") +
    # ggtitle(paste0(tmp$params$CommonName, " (", tmp$params$CombinedLatin, ")"),
    #         subtitle = "voltinism mixture model on degree-day scale") +
    labs(color = "Generation") +
    labs(x = "Degree-days accumulated (5/30°C thresholds)") +
    labs(y = "\nObserved and imputed weekly counts") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(legend.position = "none")
  mmplt
  # 
  ggsave(filename = "fig2C.tif", path = "figures", device='tiff', dpi=600)
  # saveRDS(mmplt, file = paste0("figures/species_examples/mix_", latin, ".rds"))
  dat3$region <- factor(dat3$region, levels = c("NW", "NE", "SW", "CN"))
  
  # All species mixture model figures ----
  mmplt <- ggplot(dat3 %>% filter(GenTotal>=.25), aes(x = accumdegday, y = GenTotal, color = as.factor(gen), shape = datatype)) +
    geom_point(alpha = .3) +
    scale_color_brewer(type = "qual", palette = "Dark2") +
    facet_wrap(.~region, scales = "free_y") +
    ggtitle(paste0("Voltinism mixture model for ", tmp$params$CommonName),
            subtitle = latin) +
    labs(color = "Generation", shape = "Data type") +
    labs(x = "Degree-days accumulated (5/30C thresholds)") +
    labs(y = "Observed or imputed count") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
          plot.subtitle = element_text(face = "italic"))

  mmplt
  ggsave(filename = paste(tmp$params$CommonName, "mixmod", "GDD", "png", sep = "."),
         plot = mmplt, device = "png", path = "figures/mixture_models_all", width = 10, height = 7.5, units = "in")

  
  # Abundance and phenology by generation ----
  # using GenTotal from dat3, estimate phenology/abundance for each Site/Year/Gen
  
  dat3 <- dat3 %>% 
    ungroup() %>% 
    mutate(unqID = rownames(dat3),
           genID = paste(SiteYear, gen, sep = "_"))
  
  outlist <- list()
  for (id in unique(dat3$genID)){
    tmp2 <- dat3 %>% filter(genID == id)
    out <- tmp2[1, c("genID", "SiteID", "gen", "Year", "SiteYear", "lat", "lon", "region", "RegYear", "MeanListLength", "MeanDuration", "CommonName", "YearTotal", "ObsSurvTotal", "ImpSurvTotal", "PosObsWeeks")]
    out$trap.N = TrapezoidIndex(tmp2$DOY, tmp2$GenTotal)
    out$trap.mu.doy = weighted.mean(tmp2$DOY, tmp2$GenTotal)
    out$trap.mu.gdd = weighted.mean(tmp2$accumdegday, tmp2$GenTotal)
    out$weeks.obs.gen = length(tmp2$Week[which(tmp2$datatype == "observed" & tmp2$GenTotal >= 0.5)])
    out$weeks.imp.gen = length(tmp2$Week[which(tmp2$datatype == "imputed" & tmp2$GenTotal >= 0.5)])
    out$total.obs.gen = sum(tmp2$GenTotal[which(tmp2$datatype == "observed")])
    out$total.imp.gen = sum(tmp2$GenTotal[which(tmp2$datatype == "imputed")])
    
    outlist[[length(outlist)+1]] <- out  
  }
  
  spgens <- bind_rows(outlist) %>%
    group_by(SiteID, Year) %>%
    mutate(trap.N.allgens = sum(trap.N)) %>%
    group_by(SiteID) %>%
    mutate(nyr = length(unique(Year[which(trap.N.allgens>0)])),
           YearNum = as.numeric(as.character(Year)),
           SiteIDfact = as.factor(SiteID))
  
  datlist[[sp]] <- spgens
  
} # close species loop

# mixmod results
df <- bind_rows(splist)
# population estimates by generation
df2 <- bind_rows(datlist)

saveRDS(df, "data/mixmodcomparison.rds")
saveRDS(df2, "data/genpops.rds")



# Imputation Check ----
# comparing mixmod results with and without imputation
incl_species <- mvspecies$CommonName[which(mvspecies$final_mv_analysis == 1)]

# Eastern tailed-blue has better timescale of gdd_log instead of gdd when no imputed values
# only Wild Indigo Duskwing has different best-fit # of generations (2 instead of 3) when no imputed values used (and best model is DOY)
test1 <- readRDS("data/mixmodcomparison_noimpute.rds") %>% filter(species %in% incl_species)%>% group_by(species)%>% filter(uncertainty > 0) %>% filter(uncertainty == min(uncertainty)) 
test2 <- readRDS("data/mixmodcomparison.rds") %>% filter(species %in% incl_species)%>% group_by(species)%>% filter(uncertainty > 0) %>% filter(uncertainty == min(uncertainty)) 
test1 %>% select(species, time, gen) %>% group_by(species) %>% slice(max(gen)) %>% data.frame()
test2 %>% select(species, time, gen) %>% group_by(species) %>% slice(max(gen)) %>% data.frame()

# when all fitted with best GDD model, 
# MAE of phenology peak estimate of each generation is 8.93
test1 <- readRDS("data/mixmodcomparison_noimpute.rds") %>% filter(species %in% incl_species) %>%  group_by(species) %>% filter(uncertainty > 0, time %in% c("gdd", "gdd_log")) %>% filter(uncertainty == min(uncertainty)) 
test2 <- readRDS("data/mixmodcomparison.rds") %>% filter(species %in% incl_species) %>% group_by(species) %>% filter(uncertainty > 0, time %in% c("gdd", "gdd_log")) %>% filter(uncertainty == min(uncertainty)) 
test <- left_join(test2, test1, by = c("species", "gen"))
mean(abs(test1$mean-test2$mean))

# MAE of weights (relative generation size) is 0.0129
mean(abs(test1$prop-test2$prop))

