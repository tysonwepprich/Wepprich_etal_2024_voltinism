## Header ---------------------------
## Script name: 01_monitoring_data_prep.R
## Purpose of script: Preparing monitoring data for analysis
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## This script is sourced by others

## Volunteers organized by the Ohio Lepidopterists have contributed these observations
## They ask that the data be used for research and that the group be acknowledged in publications
## https://www.ohiolepidopterists.org/
## ---

# Load packages ---- 
library(mclust)
library(lubridate)
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(viridis)

source("code/utils.R")
theme_set(theme_bw(base_size = 18)) 

# Load data ----
data <- readr::read_csv("data/data.trim.2023.csv") %>% 
  mutate(SiteID = formatC(SiteID.x, width = 3, format = "d", flag = "0"),
         SiteDate = lubridate::mdy(SiteDate))

# Two cryptic species are not distinguished in the monitoring
data$CommonName[which(data$CommonName == "Spring/Summer Azure")] <- "Azures"

# 1995 was a pilot year
data <- data %>% 
  filter(year(SiteDate) >= 1996)

# Filter out unidentified species
allspecies <- data %>% 
  filter(CommonName %in% unique(CommonName)[1:123]) %>% 
  group_by(CommonName) %>% 
  summarise(n = sum(Total)) %>% 
  arrange(n)

# Visualize surveys over time ----
surveys <- distinct(data[, c("SeqID", "SiteID", "SiteDate", "Week")])
# hist(year(surveys$SiteDate), breaks = 25)

# Which surveys are missing? Imputation check ----
# surv <- surveys %>% 
#   mutate(Year = year(SiteDate)) %>% 
#   dplyr::select(SiteID, Week, Year) %>% 
#   mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
#   distinct()
# allsurv <- surv %>% 
#   complete(SiteID, Week, Year) %>% 
#   mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
#   filter(Week <= 30) %>% 
#   filter(SiteYear %in% unique(surv$SiteYear))
# missing <- anti_join(allsurv, surv)
# 
# # we use data from sites that complete 10 or more surveys in a year
# # are missing surveys distributed differently from these sites?
# sy10cutoff <- surv %>% 
#   group_by(SiteYear) %>% 
#   mutate(nsurv = length(unique(Week))) %>% 
#   filter(nsurv >= 10)
# allsurv10 <- sy10cutoff[,c("SiteID", "Week", "Year")] %>% 
#   complete(SiteID, Week, Year) %>% 
#   mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
#   filter(Week <= 30) %>% 
#   filter(SiteYear %in% unique(sy10cutoff$SiteYear))
# missing10cutoff <- anti_join(allsurv10, sy10cutoff)
# 
# # Fig. S12 ----
# # visualize missing surveys for imputation check
# ggplot(sy10cutoff %>% filter(Week <= 30), aes(x = Week)) +
#   geom_bar(alpha = .3) +
#   # geom_bar(data = missing10cutoff, aes(x = Week), alpha = .3) +
#   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
#   ggtitle("Observed surveys by week of monitoring")
# 
# ggplot(missing10cutoff %>% filter(Week <= 30), aes(x = Week)) +
#   geom_bar(alpha = .3) +
#   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
#   ggtitle("Missing surveys by week of monitoring")


# Covariates for surveys ----
# Listlength is # of species observed, often used as catch-all covariate for effort/weather/season/etc.
covdata <- data %>%
  group_by(SeqID) %>%
  summarise(listlength = length(which(unique(CommonName) %in% allspecies$CommonName)),
            temperature = mean(c(StartTemp, EndTemp), na.rm = TRUE),
            duration = duration[1]) %>%
  distinct() %>% 
  left_join(surveys)

# impute missing duration values (due to misaligned time entries, etc.)
# only one still is NA because never reported end time for surveys.
covdata <- covdata %>% 
  mutate(Year = year(SiteDate)) %>% 
  group_by(SiteID, Year) %>% 
  mutate(duration = ifelse(is.na(duration) == TRUE | duration == 0,
                           median(duration, na.rm = TRUE),
                           duration))

# Geographic coordinates for sites ----
# address/transect start (no info in data about exact transect route)
sites <- read.csv("data/OHsites2023update.txt") %>% 
  mutate(SiteID = formatC(Name, width = 3, format = "d", flag = "0"))

# Degree days from 00_degree_days.R ----
gdd_all <- readRDS("data/daily_weather.rds")

# Summarize sites' seasons----
# # season length
# gdd_seas <- gdd_all %>% 
#   group_by(SiteID, year) %>% 
#   summarise(dd = max(accumdegday),
#             lastfdd = accumdegday[max(which(hardfrost == TRUE & yday < 170))],
#             firstfdd = accumdegday[min(which(hardfrost == TRUE & yday > 170))],
#             lastfday = yday[max(which(hardfrost == TRUE & yday < 170))],
#             firstfday = yday[min(which(hardfrost == TRUE & yday > 170))],            
#             seaslength = firstfday - lastfday,
#             seaslengthdd = firstfdd - lastfdd) %>% 
#   group_by(SiteID) %>% 
#   summarise_all(.funs = list(mean, sd), na.rm = TRUE)

# Ohio sites group into 4 regions ----
siteGDD <- gdd_all %>%
  group_by(SiteID, lat, lon) %>% 
  filter(yday == 365) %>%
  summarise(meanGDD = mean(accumdegday))

# many ways to cluster sites, but using lat/lon is simplest 
# wanted 4 regions for plotting simplicity
sitemod <- densityMclust(siteGDD[,c(2,3)], G = 1:4, modelNames = "EVV", plot = FALSE)
siteGDD$region <- as.character(sitemod$classification)

# # visualize regional clusters, imagine an Ohio outline
# a <- ggplot(data = siteGDD, aes(x = lon, y = lat, group = region, color = region)) + 
#   geom_point(size = 2)
# a

siteGDD$region <- plyr::mapvalues(siteGDD$region, from = c("1", "2", "3", "4"), 
                                  to = c("NE", "NW", "CN", "SW"))
gdd <- gdd_all %>% 
  left_join(siteGDD[, c("SiteID", "region")]) %>% 
  mutate(SiteID = formatC(SiteID, width = 3, format = "d", flag = "0"),
         SiteDate = paste(year, yday, sep = "-"),
         SiteDate = as.Date(parse_date_time(SiteDate, orders = "Yj")))

# Mean site # yrs and # surveys/yr ----
# # Reported in Methods
# samplesize <- surveys %>% mutate(year = year(SiteDate)) %>% group_by(SiteID, year) %>% 
#   summarise(nsurv = length(unique(SeqID))) %>% 
#   ungroup() %>% 
#   group_by(SiteID) %>% 
#   summarise(nyr = length(unique(year)),
#             meansurv = mean(nsurv))
# summary(samplesize)
