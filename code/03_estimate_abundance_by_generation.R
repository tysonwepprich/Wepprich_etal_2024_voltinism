## Header ---------------------------
## Script name: 03_estimate_abundance_by_generation.R
## Purpose of script: Uses GAM to impute missing surveys, fits mixture models, and estimates abundance/phenology for each generation
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## 
## ---

source('01_data_prep.R')

# 42 candidate species for mixture models, not all will work (mv_analysis == 1)
# After species filter (1st pass with mixture models) 35 remain (final_mv_analysis == 1)
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
  minvolt <- mvspecies$MinVoltMM[sp]
  maxvolt <- mvspecies$Voltinism[sp]
  
  # import GAM phenology model and count data used for it
  mf <- modfiles[grep(pattern = species, x = modfiles, fixed = TRUE)]
  tmp <- readRDS(mf)
  mod <- tmp$gammod
  dat <- tmp$datGAM
  
  # find missing weeks to impute
  
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
  
  # Predict missing weeks. Better to round to whole number or simulate from GAM's negbinom distribution?
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
  
  # Mclust docs: what's the uncertainty measure? Could be good way to compare DOY/GDD and regional split
  dd_dist <- rep(alldat$accumdegday, round(alldat$Total))
  dd_doy <- rep(alldat$DOY, round(alldat$Total))
  dd_log <- rep(alldat$accumdegday, round(log(alldat$Total + 1)))
  
  mm <- Mclust(dd_dist, G = minvolt:maxvolt, modelNames = "E")
  mm_doy <- Mclust(dd_doy, G = minvolt:maxvolt, modelNames = "E")
  mm_log <- Mclust(dd_log, G = minvolt:maxvolt, modelNames = "E")
  # mm2 <- Mclust(cbind(dd_dist, dd_doy), G = minvolt:maxvolt)
  
 
  
  # test doy vs gdd
  # constraining to equal variance here, but could try unequal for more realism
  mmdf <- data.frame(species = species, gen = 1:length(mm$parameters$mean),
                     prop = mm$parameters$pro, mean = mm$parameters$mean,
                     time = "gdd", uncertainty = mean(mm$uncertainty))
  mmdf2 <- data.frame(species = species, gen = 1:length(mm_doy$parameters$mean),
                      prop = mm_doy$parameters$pro, mean = mm_doy$parameters$mean,
                      time = "doy", uncertainty = mean(mm_doy$uncertainty))
  mmdf3 <- data.frame(species = species, gen = 1:length(mm_log$parameters$mean),
                      prop = mm_log$parameters$pro, mean = mm_log$parameters$mean,
                      time = "gdd_log", uncertainty = mean(mm_log$uncertainty))
  splist[[sp]] <- bind_rows(mmdf, mmdf2, mmdf3)
  
  # CHECK THAT PREDICTIONS MAKE SENSE
  if(mmdf$uncertainty[1]<mmdf3$uncertainty[1] & length(mm$parameters$mean)>1){
    mm <- Mclust(dd_dist, G = minvolt:maxvolt, modelNames = "E")
  }else{
    mm <- Mclust(dd_log, G = minvolt:maxvolt, modelNames = "E")
  }
  
 # assigning generation
  # avoiding rounding so that GAM predictions multiplied by mixed model generation probabilities
  # each count split same way each time
  dat <- alldat %>% 
    group_by(SiteYear) %>% 
    mutate(YearTotal = sum(Total[which(datatype == "observed")]),
           ObsSurvTotal = length(Week[which(datatype == "observed")]),
           ImpSurvTotal = length(Week[which(datatype == "imputed")]),
           PosObsWeeks = length(Week[which(datatype == "observed" & Total > 0)]),
           MeanListLength = mean(listlength, na.rm = TRUE),
           MeanDuration = mean(duration, na.rm = TRUE))
  
  mmpred <- predict(mm, newdata = dat$accumdegday)
  res <- as.data.frame(mmpred$z)
  
  dat2 <- cbind(dat, dat$Total * res)
  dat3 <- dat2 %>% 
    pivot_longer(cols = as.character(c(1:maxvolt)),
                 names_to = "gen",
                 values_to = "GenTotal")
  
  dat3$GenTotal <- round(dat3$GenTotal, digits = 3)
  
  ###################
  
  
  dat3$region <- factor(dat3$region, levels = c("NW", "NE", "SW", "CN"))
  
  mmplt <- ggplot(dat3 %>% filter(GenTotal>0), aes(x = accumdegday, y = GenTotal, color = as.factor(gen), shape = datatype)) +
    geom_point(alpha = .2) +
    scale_color_brewer(type = "qual", palette = "Dark2") +
    facet_wrap(.~region, scales = "free_y") +
    ggtitle(paste0(tmp$params$CommonName, " (", tmp$params$CombinedLatin, ")"),
            subtitle = "voltinism mixture model on degree-day scale") +
    labs(color = "Generation") +
    labs(x = "Degree-days accumulated (5/30C thresholds)") +
    labs(y = "Observed or imputed count")
  # mmplt
  ggsave(filename = paste(tmp$params$CommonName, "mixmod", "GDD", "png", sep = "."),
         plot = mmplt, device = "png", path = "gams/plots", width = 8, height = 6, units = "in")
  
  
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

#mixmod results
df <- bind_rows(splist)
#population estimates by generation
df2 <- bind_rows(datlist)

saveRDS(df, "mixmodcomparison.rds")
saveRDS(df2, "genpops.rds")


# # added ETB
# df <- readRDS("mixmodcomparison.rds")
# df <- bind_rows(df, splist[[13]])
# saveRDS(df, "mixmodcomparison.rds")
# df2 <- readRDS("genpops.rds")
# df2 <- bind_rows(df2, datlist[[13]])
# saveRDS(df2, "genpops.rds")

# this code was to visualize GDD vs DOY comparison
# df2 <- df %>% 
#   group_by(species) %>%
#   arrange(time) %>% 
#   summarise(maxgen = max(gen),
#             unratio = uncertainty[1] / uncertainty[length(uncertainty)])
# 
# df2 <- df %>% 
#   pivot_wider(id_cols = c(species, gen), names_from = time, values_from = prop)
# 
# ggplot(df2, aes(x = gdd, y = doy, color = as.factor(gen))) +
#   geom_point(size = 3) +
#   geom_abline(slope = 1, intercept = 0) +
#   facet_wrap(~species)


# Remove some species from analysis, looking at GAMs and mix models
# App Brown (uni), 
# Azures (multispecies), 
# Clouded Sulphur (end of season),
# Dun Skipper (iftner says 3 broods?), 
# ETB (too blended)
# ETS (3rd maybe not there, <10% in GDD mixmid), 
# Gray HS (3-4 iftner, all mixmod miss early 1st generation),
# Monarch (but try on the side), 

# Check
# Common Sooty (only 2 in DOY mixmod, 3 in GDD mixmod & GAM)
# Euro Skip (2nd maybe not there, <7% in GDD Log mixmod), 
# Silvery Check (should have 2 not 3 gens, adjust maxgen?),
# WIDW (should be 3, GDD mixmod selected 2 bc missed small 1st gen, LOG mixmod correct, GAM obviously 3)

# 2nd filter problems, solved by using gdd_log mixmod!
# Common Sooty wrong peaks selected in mixmod
# Euroskip double 1st gen selected
# Viceroy should be 3 gen (2 selected)
# WIDW (catches 1st gen if regionally stratified)

