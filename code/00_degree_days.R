## Header ---------------------------
## Script name: 00_degree_days.R
## Purpose of script: Download weather data for each site from Daymet
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## Calculate daily degree-days with 5C/30C thresholds
## Calculate photoperiod for each day
## Calculate days with hard frosts (<=-2C)
## Weather data from https://daymet.ornl.gov/
## ---

library(dplyr)
library(daymetr)
source("code/utils.R") # for degree-day and photoperiod functions

# Long-term monitoring sites from Ohio Lepidopterists
sites <- read.csv("data/OHsites2023update.txt") %>% 
  mutate(SiteID = formatC(Name, width = 3, format = "d", flag = "0"))

startyear <- 1995
endyear <- 2022

# for dataframe of sites
siteslist <- list()
for(i in 1:nrow(sites)){
  temp <- download_daymet(site = sites$SiteID[i], lat = sites$lat[i], lon = sites$lon[i],
                          start = startyear, end = endyear, internal = TRUE,
                          silent = TRUE, force = FALSE)
  outdf <- temp$data %>%
    mutate(elev = temp$altitude,
           SiteID = sites$SiteID[i],
           lat = sites$lat[i],
           lon = sites$lon[i])
  siteslist[[i]] <- outdf
}

gdd <- bind_rows(siteslist)

# lower developmental threshold of 5 instead of 10 so that April surveys differentiated
ldt <- 5
udt <- 30

gdd_all <- gdd %>%
  mutate(SiteID = formatC(SiteID, width = 3, format = "d", flag = "0")) %>%
  rowwise() %>%
  mutate(degday = degree_days(tmin..deg.c., tmax..deg.c., ldt, udt, method = "single.sine")) %>%
  ungroup() %>%
  group_by(SiteID, year) %>%
  arrange(yday) %>%
  mutate(SiteDate = as.POSIXct(strptime(paste(year, yday, sep = " "), format = "%Y %j")),
         accumdegday = cumsum(degday),
         daylength = photoperiod(lat, yday),
         lightfrost = tmin..deg.c. <= 0,
         hardfrost = tmin..deg.c. <= -2) %>% 
  dplyr::select(year, yday, tmax..deg.c., tmin..deg.c., SiteID:hardfrost)

saveRDS(gdd_all, "data/daily_weather.rds")


