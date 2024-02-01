## Header ---------------------------
## Script name: 02_fit_plot_gam.R
## Purpose of script: Fit seasonal phenological models using generalized additive models
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## GAM formula is same as Wepprich et al. (2019) https://doi.org/10.1371/journal.pone.0216270
## Output is too big to save in repository (each species file up to 900 MB)
## Fitting models for 30 species takes 45 minutes on a 8-core 32GB-RAM laptop
## Running this script is not necessary to reproduce the main results, but would
## be required to run 03_estimate_abundance_by_generation.R to fit mixture models with GAM-imputed counts
## ---

source('code/01_data_prep.R')

# Fit a generalized additive model (GAM) for each species' counts
# This model will be used to impute missing counts a la UKBMS
# Also will plot seasonal phenological patterns in different regions/years

# Select species, as is it will attempt to fit models for 100+ species

# run examples species only (Phyciodes tharos) for Fig 2 or all 30 species
allspecies <- allspecies[which(allspecies$CommonName == "Pearl Crescent"), ]
# allspecies <- read.csv("data/speciestraits.csv") %>% 
  # filter(final_mv_analysis == 1) %>% 
  # dplyr::select(CommonName, final_mv_analysis)


system.time({
# Fit GAM ----
# use parallel processing to fit GAMs for each species
# results in large files saved
library(parallel)
library(foreach)
ncores <- 8
if(ncores > (parallel::detectCores() / 2)){
  ncores <- parallel::detectCores() / 2
}
cl <- makePSOCKcluster(ncores)
doParallel::registerDoParallel(cl)

mcoptions <- list(preschedule = FALSE)

# foreach loop
outfiles <- foreach(sp = 1:nrow(allspecies),
                    .combine='c',
                    .packages= c("mgcv", "dplyr", "tidyr", "purrr",
                                 "lubridate", "mclust"),
                    .export = c("data", "surveys", "covdata", "allspecies", "gdd"),
                    .inorder = FALSE,
                    .options.multicore = mcoptions) %dopar% {
                      
                      species <- allspecies$CommonName[sp]
                      pars <- allspecies %>% filter(CommonName == species)
                      
                      counts <- data %>% 
                        filter(CommonName == species) %>% 
                        mutate(DOY = yday(SiteDate),
                               Year = year(SiteDate))
                      
                      #get unique surveys, including those where species not counted (but were at one time)
                      survs <- surveys %>% 
                        # filter(year(SiteDate) %in% unique(counts$Year)) %>% 
                        filter(SiteID %in% unique(counts$SiteID)) %>% 
                        mutate(Year = year(SiteDate))
                      
                      #Add zeros to surveys when species not counted during a survey
                      test <- left_join(survs, counts, by = c("SiteID", "SiteDate", "Week", "SeqID", "Year"))
                      counts <- test[, c("SiteID", "SiteDate", "Week", "SeqID", "Total", "Year")]
                      counts$Total[which(is.na(counts$Total))] <- 0 
                      counts <- left_join(counts, covdata)
                      counts$temperature[which(counts$temperature < 50)] <- NA
                      counts$duration[which(counts$duration == 0)] <- NA
                      
                      # scaling covariates
                      # previously, duration and temperature not that important
                      # will just use list-length scaled by site and month
                      # this accounts for survey variation while controlling for site and season
                      counts <- counts %>% 
                        group_by(SiteID) %>% 
                        mutate(zlistlength = as.numeric(scale(listlength)))
                      
                      
                      # trying to add GDD instead of ordinal date
                      counts <- counts %>% 
                        left_join(gdd) %>% 
                        group_by(SiteID, Year) %>% 
                        mutate(SurvPerYear = length(unique(SeqID)),
                               YearTotal = sum(Total))
                      
                      # what if cutoff for inclusion is really open?
                      # could filter out sites with lower effort later
                      # also this better accounts for SiteYear zero counts, 
                      # which are important for population trends
                      dat <- counts %>% filter(YearTotal >= 0, SurvPerYear >= 5)
                      
                      mod <- list()
                      
                      if(nrow(dat) == 0){
                        mod$error <- "no data"
                      }else{
                        
                        dat$Year <- as.factor(as.character(dat$Year))
                        dat$region <- as.factor(as.character(dat$region))
                        dat$SiteID <- as.factor(as.character(dat$SiteID))
                        dat$SiteYear <- as.factor(paste(dat$SiteID, dat$Year, sep = "_"))
                        dat$zlistlength[which(is.na(dat$zlistlength))] <- 0
                        dat$RegYear <- as.factor(paste(dat$region, dat$Year, sep = "_"))
                        dat$DOY <- dat$yday
                        dat <- as.data.frame(dat)
                        
                        dat <- dat[which(!is.na(dat$accumdegday)), ]
                        
                        
                        if(sum(dat$Total) < 20|length(unique(dat$SiteID)) < 2|length(unique(dat$Year)) < 2|
                           length(unique(dat$SiteYear))<2|length(unique(dat$RegYear))<2) {
                          mod$error <- "few data"
                        }else{
                          safe_bam <- purrr::safely(bam)
                          
                          modtime_nb <- system.time({ 
                            mod_nb <- safe_bam(Total ~
                                                 # s(zlistlength, bs = "cr") +
                                                 te(lat, lon, accumdegday, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                                                 s(SiteYear, bs = "re") +
                                                 s(Year, bs = "re") +
                                                 s(SiteID, bs = "re") +
                                                 s(RegYear, DOY, bs = "fs", k = 5, m = 1),
                                               family = nb(theta = NULL, link = "log"),
                                               # family = poisson(link = "log"),
                                               data = dat,
                                               method = "fREML",
                                               discrete = TRUE, nthreads = 2,
                                               control = list(maxit = 500))
                            
                          })
                          
                          modtime_po <- system.time({ 
                            mod_po <- safe_bam(Total ~
                                                 # s(zlistlength, bs = "cr") +
                                                 te(lat, lon, accumdegday, bs = c("tp", "cr"), k = c(5, 30), d = c(2, 1)) +
                                                 s(SiteYear, bs = "re") +
                                                 s(Year, bs = "re") +
                                                 s(SiteID, bs = "re") +
                                                 s(RegYear, DOY, bs = "fs", k = 5, m = 1),
                                               # family = nb(theta = NULL, link = "log"),
                                               family = poisson(link = "log"),
                                               data = dat, 
                                               method = "fREML",
                                               discrete = TRUE, nthreads = 2,
                                               control = list(maxit = 500))
                            
                          })
                          
                          if(is.null(mod_po$error) == TRUE & is.null(mod_nb$error) == TRUE){
                            # compare models
                            if(AIC(mod_po$result) <= AIC(mod_nb$result)){
                              mod <- mod_po
                              modtime <- modtime_po
                            }else{
                              mod <- mod_nb
                              modtime <- modtime_nb  
                            }
                          }else{
                            if(is.null(mod_po$error) == TRUE){
                              mod <- mod_po
                              modtime <- modtime_po
                            }else{
                              mod <- mod_nb
                              modtime <- modtime_nb
                            }
                          }
                        }
                      }
                      
                      if(is.null(mod$error)){
                        pars$modtime <- as.numeric(modtime)[1]
                        pars$AIC <- AIC(mod$result)
                        summod <- summary(mod$result)
                        pars$N <- summod$n
                        pars$dev.expl <- summod$dev.expl
                        pars$family <- summod$family$family
                        outlist <- list()
                        outlist[["params"]] <- pars
                        outlist[["gammod"]] <- mod$result
                        outlist[["gamtime"]] <- as.numeric(modtime)[1]
                        outlist[["datGAM"]] <- dat
                        saveRDS(outlist, paste(paste0("gams/", species), "5", "final", "rds", sep = "."))
                        
                      }else{
                        outlist <- list()
                        outlist[["params"]] <- pars
                        outlist[["gammod"]] <- mod$error
                        outlist[["datGAM"]] <- dat
                        saveRDS(outlist, paste(paste0("gams/", "gamerr"), species, "5", "rds", sep = "."))
                      }
                      return(sp)
                    }

if(.Platform$OS.type == "windows"){
  stopCluster(cl)
}
})

# Plot GAM ----
# Plot seasonal phenological variation from species GAM
## ---
## Notes: 
## Running this script is not necessary to reproduce the main results
## Uses output from 02_fit_gam.R that is not included in repository due to file size
## These plots were used to identify candidate multivoltine species
## Included in example plot of Pearl Crescent (Fig 2)

source('code/01_data_prep.R')

modfiles <- list.files(path = "gams", full.names = TRUE)
modfiles <- modfiles[grep(pattern = "final", x = modfiles, fixed = TRUE)]
sp_params <- list()

for (mf in seq_along(modfiles)){
  
  tmp <- readRDS(modfiles[mf])
  mod <- tmp$gammod
  counts <- tmp$datGAM
  sp_params[[mf]] <- tmp$params
  
  # name for figures
  latin <- read.csv("data/species_names.csv") %>% 
    dplyr::select(CommonName, Genus, Species) %>%
    filter(CommonName == tmp$params$CommonName) %>% 
    mutate(Latin = paste(Genus, Species, sep = " ")) %>% 
    pull(Latin) %>% 
    as.character()
  
  # plot GAM predictions for species/model
  preds <- gdd %>% 
    mutate(Year = year) %>% 
    mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
    filter(SiteYear %in% unique(counts$SiteYear)) %>%
    group_by(SiteID, Year) %>% 
    mutate(SiteYearGDD = max(accumdegday),
           DOY = yday) %>% 
    filter(DOY %in% seq(90, 305, 1)) %>% 
    ungroup() %>% 
    mutate(RegYear = as.factor(paste(region, Year, sep = "_" )),
           Year = as.factor(as.character(Year)),
           SiteID = as.factor(SiteID),
           SiteYear = as.factor(SiteYear))
  
  preds$adjY <- predict.gam(object = mod, newdata = preds, type="response")
  
  preds2 <- preds %>% 
    group_by(SiteYear) %>%
    mutate(Gamma = adjY / as.vector(sum(adjY)),
           SiteYearTotal = sum(adjY)) %>%
    ungroup() %>% 
    filter(adjY > 0, year %in% c(2008:2017), SiteYearTotal > median(SiteYearTotal)) %>% 
    mutate(yrgroup = ifelse(year %in% c(2010, 2011, 2012, 2016, 2017), "Warm years", "Cool years"),
           sitegroup = ifelse(region %in% c("NE", "NW"), "Cool sites", "Warm sites")) %>% 
    group_by(sitegroup, yrgroup) %>% 
    filter(SiteYear %in% sample(unique(SiteYear), 5, replace = TRUE))
  
  # preds$region <- factor(preds$region, levels = c("NW", "NE", "SW", "CN"))
  
  # Fig 2A ----
  # Degree-day timescale
  gamplt <- ggplot(preds2, aes(x = accumdegday, y = Gamma, group = SiteYear, color = yrgroup)) +
    geom_path(alpha = .75) + 
    scale_colour_brewer(name = NULL, palette = "Set1", direction = -1) +   
    facet_wrap(~sitegroup, ncol = 1, scales = "free_y") +
    # ggtitle(paste0(tmp$params$CommonName, " (", latin, ")"),
            # subtitle = "Seasonal phenology modeled on degree-day scale") +
    labs(x = "Degree-days accumulated (5/30°C thresholds)") +
    labs(y = "Scaled phenology (model predictions)") +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(legend.position = c(.15, .925))
  gamplt
  ggsave(filename = "example_gam.tif", path = "figures", device='tiff', dpi=600)
  
  saveRDS(gamplt, file = paste0("figures/species_examples/gam_", latin, ".rds"))
  # ggsave(filename = paste(tmp$params$CommonName, "GAM", "GDD", "png", sep = "."), 
         # plot = gamplt, device = "png", path = "figures/species_examples/", width = 8, height = 6, units = "in")
  
  # Fig 2B ----
  # Day of year timescale
  preds2$date <- as.Date(as.Date("2019-12-31",  "%Y-%m-%d") + ddays(preds2$yday))
  
  gamdoyplt <- ggplot(preds2, aes(x = date, y = Gamma, group = SiteYear, color = yrgroup)) +
    geom_path(alpha = .75) + 
    scale_colour_brewer(name = NULL, palette = "Set1", direction = -1) +   
    scale_x_date(date_breaks = "1 month", date_labels = "%b") +
    facet_wrap(~sitegroup, ncol = 1, scales = "free_y") +
    # ggtitle(paste0(tmp$params$CommonName, " (", latin, ")"),
    # subtitle = "Seasonal phenology modeled on degree-day scale") +
    labs(x = "Calendar date") +
    labs(y = "Scaled phenology (model predictions)") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
    theme(legend.position = c(.15, .925))
  gamdoyplt
  ggsave(filename = "example_gamdoy.tif", path = "figures", device='tiff', dpi=600)
  
  
  
}

outparams <- bind_rows(sp_params)
saveRDS(outparams, "data/modparams.5.rds")


# GAM parameters ----
# This collects information from the large output from the GAMs
# Used to report on model fit (deviance explained) in Results
# need "sp" to limit results to species in analysis
# "sp" is a vector of species from 07_species_models.R once final list for inclusion was determined.

modfiles <- list.files(path = "gams", full.names = TRUE)
modfiles <- modfiles[grep(pattern = "final", x = modfiles, fixed = TRUE)]
sp_params <- list()

for (mf in seq_along(modfiles)){
  
  tmp <- readRDS(modfiles[mf])
  sp_params[[mf]] <- tmp$params
}
outparams <- bind_rows(sp_params)

outparams <- readRDS("data/modparams.5.rds") %>% 
  filter(CommonName %in% sp)





