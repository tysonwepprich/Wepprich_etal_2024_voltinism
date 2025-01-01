## Header ---------------------------
## Script name: 04_prep_model_data.R
## Purpose of script: Creating variables for last generation size and weather covariates
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## Running this script is not necessary to reproduce the main results, because 
## output saved "data/modeling_data.rds" for use in later scripts
## Also includes simple models of weather variables over time (reported in results).
## Some visualizations of weather and last generation variation (not included in paper)
## 
## ---

source('code/01_data_prep.R')
library(broom.mixed)
library(lme4)

mvspecies <- read.csv("data/speciestraits.csv", header = TRUE) %>%
  # filter(mv_analysis == 1) %>%
  filter(final_mv_analysis == 1) %>%
  droplevels.data.frame()
incl_species <- mvspecies$CommonName[which(mvspecies$final_mv_analysis == 1)]


genpop <- readRDS("data/genpops.rds") %>% 
  filter(CommonName %in% incl_species) 

# Filter data, add last generation size variables ----
popdat <- genpop %>% 
  filter(nyr >= 5, ObsSurvTotal >= 10, PosObsWeeks >= 3, trap.N.allgens >= 10, YearTotal >= 3) %>% 
  mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
  group_by(CommonName) %>% 
  mutate(maxbrood = base::max(as.numeric(gen))) %>% 
  group_by(CommonName, SiteYear) %>%
  arrange(gen) %>%
  mutate(lastprop = trap.N[maxbrood] / (trap.N[maxbrood-1] + trap.N[maxbrood]),
         lastN = trap.N[maxbrood],
         lastDOY = trap.mu.doy[maxbrood],
         lastGDD = trap.mu.gdd[maxbrood],
         lastratio = trap.N[maxbrood] / trap.N[maxbrood-1],
         lastloglambda = log((trap.N[maxbrood]+1) / (trap.N[maxbrood-1]+1)),
         lastObs = total.obs.gen[maxbrood],
         lastImp = total.imp.gen[maxbrood])


# Compare quantification of voltinism
# natural log of ratio of last / penultimate works for models of population growth rates
# main area of mismatch is when last generation is very small. 
# lastprop would be near 0 no matter how large the penultimate generation, but
# with lastloglambda a last generation of near 0 could have varying "sizes" depending 
# on growth rate from penultimate generation size.
ggplot(popdat, aes(x = lastprop, y = lastloglambda, color = log(lastN+1))) +
  geom_point() +
  facet_wrap(~CommonName)
ggplot(popdat, aes(x = lastloglambda)) +
  geom_density() +
  facet_wrap(~CommonName)



# Site covariates ----
sitecoords <- sites[,c("lon", "lat")] 

gdd$SiteYear <- paste(gdd$SiteID, gdd$year, sep = "_")

# what about expected gdd by photoperiod/site?
# gddmismatch > 0 means that there are more gdd left in the season than expected by the site/date average
gdd_left <- gdd %>% 
  group_by(SiteID, year) %>% 
  mutate(siteyrtotalgdd = max(accumdegday),
         firstfrost = min(accumdegday[which(yday > 200 & hardfrost == TRUE)])) %>% 
  ungroup() %>% 
  mutate(actualgddleft = siteyrtotalgdd - accumdegday,
         actualgddleftfrost = firstfrost - accumdegday) %>% 
  group_by(SiteID, yday) %>% 
  mutate(expgddleft = mean(actualgddleft),
         expgddleftfrost = mean(actualgddleftfrost)) %>% 
  ungroup() %>% 
  mutate(gddmismatch = actualgddleft - expgddleft,
         gddmismatchfrost = actualgddleftfrost - expgddleftfrost) 


# Covariates at brood peaks ----
# UNCOMMENT FOR FIRST RUN
# 
# # Run once, took 15 minutes
# # get photoperiod at brood peaks, gdd left in year at peak, temperature around peak
# gdd2vars <- function(SiteYr, meanmu){
#   temp <- gdd_left %>%
#     filter(SiteYear == SiteYr) %>%
#     mutate(daymeantemp = (tmax..deg.c. + tmin..deg.c.) / 2)
# 
# 
#   ord     <-     temp$yday[which(temp$accumdegday > meanmu)[1]]
#   photo   <-     temp$daylength[which(temp$accumdegday > meanmu)[1]]
#   gddleft <-     temp$actualgddleft[which(temp$accumdegday > meanmu)[1]]
#   gddexpected <- temp$expgddleft[which(temp$accumdegday > meanmu)[1]]
#   gddmismatch <- temp$gddmismatch[which(temp$accumdegday > meanmu)[1]]
#   gddleftfrost <-     temp$actualgddleftfrost[which(temp$accumdegday > meanmu)[1]]
#   gddexpectedfrost <- temp$expgddleftfrost[which(temp$accumdegday > meanmu)[1]]
#   gddmismatchfrost <- temp$gddmismatchfrost[which(temp$accumdegday > meanmu)[1]]
# 
#   return(data.frame(ord,
#              photo,
#              gddleft,
#              gddexpected,
#              gddmismatch,
#              gddleftfrost,
#              gddexpectedfrost,
#              gddmismatchfrost))
# }
# 
# 
# system.time({covs <- popdat %>%
#   group_by(CommonName, SiteYear, gen) %>%
#   do(gdd2vars(.$SiteYear, .$trap.mu.gdd))
# })
# saveRDS(covs, "data/covs.rds")

covs <- readRDS("data/covs.rds")

alldat <- popdat %>% 
  left_join(covs)

# modeling data will contain penultimate generation, last generations, and following 1st generation
# penultimate
moddat <- alldat %>% 
  group_by(CommonName) %>%
  mutate(uniqSY = length(unique(SiteYear))) %>% 
  filter(gen == maxbrood-1, trap.N > 1)

# add following 1st generation data for model of overwinter population growth rate
# size/time of penult gen, size/time of last gen and size/time of 1st gen next year
firstdat <- alldat %>% 
  filter(gen == 1) %>% 
  group_by(CommonName, SiteID) %>% 
  arrange(YearNum) %>% 
  mutate(firstN = lead(trap.N, 1),
         firstDOY = lead(trap.mu.doy, 1),
         firstGDD = lead(trap.mu.gdd, 1),
         firstYear = lead(YearNum, 1),
         firstObs = lead(total.obs.gen, 1),
         firstImp = lead(total.imp.gen, 1)) %>% 
  filter(firstYear == (YearNum + 1)) %>% 
  filter(firstN > 1) %>% 
  dplyr::select(SiteID, Year, SiteYear, CommonName, firstN:firstImp)


moddat <- moddat %>% 
  left_join(firstdat)

# Winter covariates ----
# Used 5 month winter mean minimum because it matches the months without butterfly monitoring
winter <- gdd %>% 
  group_by(SiteID) %>% 
  arrange(SiteDate) %>% 
  mutate(mean3month = zoo::rollmean(tmin..deg.c., 90, fill = NA, align = "left"),
         mean5month = zoo::rollmean(tmin..deg.c., 150, fill = NA, align = "left")) %>% 
  group_by(SiteYear, year, SiteID) %>% 
  summarise(falltemp = mean3month[which(yday == 244)],
            wintertemp = mean3month[which(yday == 365)],
            allwintertemp = mean5month[which(yday == 305)],
            siteyrtotalgdd = max(accumdegday),
            firstfrostgdd = min(accumdegday[which(yday > 200 & hardfrost == TRUE)]),
            firstfrostdoy = min(yday[which(yday > 200 & hardfrost == TRUE)]))

moddat <- moddat %>% 
  left_join(winter)

# filter out species with few available lambdas
moddat <- moddat %>% 
  group_by(CommonName) %>% 
  mutate(uniqSYlam = length(unique(SiteYear[which(!is.na(firstN))])))

moddat <- moddat %>% 
  filter(uniqSY >= 40 & uniqSYlam > 20) # filter from 35 to 30 species

saveRDS(moddat, "data/modeling_data.rds")

# Fig 2D ----
# # Original
# # example species
# lgvar <- moddat %>%
#   filter(CommonName == "Pearl Crescent", year %in% c(2008:2017), YearTotal >= 50, nyr >= 15, photo > 12.2) %>%
#   mutate(yrgroup = ifelse(year %in% c(2010, 2011, 2012, 2016, 2017), "Warm years", "Cool years"),
#          sitegroup = ifelse(region %in% c("NE", "NW"), "Cool sites", "Warm sites"))
# lgvarplt <- ggplot(data = lgvar, aes(x = photo,
#                                       y = lastratio,
#                                       color = yrgroup))+
#   geom_point()+
#   scale_x_reverse() +
#   scale_y_log10() +
#   scale_color_brewer(name = NULL, palette = "Set1", direction = -1) +
#   scale_fill_brewer(name = NULL, palette = "Set1", direction = -1) +
#   facet_grid(sitegroup~yrgroup) +
#   theme_bw(base_size = 18) +
#   theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
#         legend.position = c(.25, .1)) +
#   labs(
#     x = "Day length at penultimate peak",
#     y = "Last generation size",
#     title = NULL)
# lgvarplt
# 
# ggsave(filename = "Fig2D.tif", path = "figures", device='tiff', dpi=600)

# revision Fig. 2D ----
lgvar <- moddat %>% 
  filter(CommonName == "Pearl Crescent", 
         year %in% c(2008:2017),
         YearTotal >= 50, nyr >= 15, photo > 12.2) %>%
  mutate(yrgroup = ifelse(year %in% c(2010, 2011, 2012, 2016, 2017), "Warm years", "Cool years"),
         sitegroup = ifelse(region %in% c("NE", "NW"), "Cool sites", "Warm sites")) %>% 
  group_by(SiteID) %>% 
  mutate(site_photo = mean(photo),
         site_LG = DescTools::Gmean(lastratio),
         site_gdd = mean(siteyrtotalgdd))
lgsite <- lgvar %>% select(SiteID, site_photo, site_LG) %>% distinct()


lgvarplt <- ggplot(data = lgvar, aes(x = photo, y = lastratio, group = SiteID, color = site_gdd)) +
  # geom_point(alpha = .1) +
  geom_line(stat="smooth", method = "lm", se = FALSE,
            size = 1,
            # linetype ="dashed",
            alpha = 0.5) +  
  scale_color_gradient(name = "Mean site\ndegree-days", high = "#ca0020", low = "#0571b0", breaks = c(2710, 3380)) +
  scale_y_log10() +
  scale_x_reverse() +
  geom_point(data = lgsite, aes(x = site_photo, y = site_LG), inherit.aes = FALSE, color = "black", size = 1.5) +
  geom_smooth(data = lgsite, aes(x = site_photo, y = site_LG), inherit.aes = FALSE, method = "lm", se = FALSE, color = "black", size = 1.25) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = c(.85, .75)) +
  labs(
    x = "Day length (hours) at mean phenology\nof penultimate generation",
    y = "Relative size of last generation compared\nto penultimate generation size (log scale)",
    title = NULL)

lgvarplt
ggsave(filename = "Fig2D.tif", path = "figures", device='tiff', dpi=600)


# # good figure of raw data
# # clear to see some trends (and outliers from misclassifying broods)
# # illustrates that some species vary more in voltinism across regions or across years
# summdat <- moddat %>% 
#   group_by(CommonName, Year, region) %>% 
#   summarise(size = log(sum(trap.N)),
#             mu.gdd=mean(trap.mu.gdd),
#             latitude=mean(lat),
#             ord=mean(ord),
#             nsites=length(unique(SiteID)),
#             photo=mean(photo),
#             temp=mean(meantemp),
#             gddmm=mean(gddmismatch),
#             lastprop=mean(lastprop),
#             lastlam = mean(lastloglambda))
# 
# 
# vistest <- ggplot(data = summdat, aes(x = ord, 
#                                       y = lastlam, 
#                                       group = CommonName, 
#                                       color = region))+
#   geom_point()+
#   facet_wrap(~CommonName, scales = "free")+
#   theme_bw()
# vistest


# Environmental covariates ----
# Correlations among environmental covariates

# Example using Pieris rapae
# using degree-days left in the season isn't useful because it's so correlated with the last generation size
envdat <- moddat %>% 
  filter(CommonName == "Cabbage White")
GGally::ggpairs(envdat[, c("lastprop", "gddleft", "gddexpected", "gddleftfrost", "gddmismatch", "allwintertemp")])
# GDD left, expected correlated highly with lastprop. Makes sense, all caused by earlier penult phenology.
GGally::ggpairs(envdat[, c("lastprop", "gddmismatch","gddmismatchfrost", "firstfrostdoy", "firstfrostgdd", "allwintertemp")])
# firstfrostdoy are not correlated with other variables like lastprop or allwintertemp 
GGally::ggpairs(gdd_left[gdd_left$yday == 225, c("gddmismatch","gddmismatchfrost", "firstfrost", "actualgddleft")])

# frost is more stochastic than gddmismatch (more variation across sites within a particular year)
ggplot(envdat, aes(x = YearNum, y = firstfrostdoy)) +
  geom_jitter() + geom_smooth() # + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal()
ggplot(envdat, aes(x = YearNum, y = gddmismatch)) +
  geom_jitter() + geom_smooth() # + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal()

# GDD trends comparing before last gen peak and after last gen peak
# TODO: models to quantify these trends, SiteID intercepts needed
ggplot(moddat, aes(x = YearNum, y = trap.mu.gdd)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~CommonName)
ggplot(moddat, aes(x = YearNum, y = trap.mu.gdd, group = CommonName)) +
  geom_smooth(method = "lm")
ggplot(moddat, aes(x = YearNum, y = gddleftfrost)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~CommonName)
ggplot(moddat, aes(x = YearNum, y = gddleftfrost, group = CommonName)) +
  geom_smooth(method = "lm")

# Weather trends ----
# Reported in results: environmental covariates

weather_trends <- gdd %>% 
  group_by(SiteID) %>% 
  arrange(SiteDate) %>% 
  mutate(mean3month = zoo::rollmean((tmin..deg.c. + tmax..deg.c.)/2, 90, fill = NA, align = "left"),
         mean5month = zoo::rollmean((tmin..deg.c. + tmax..deg.c.)/2, 150, fill = NA, align = "left")) %>% 
  group_by(SiteYear, year, SiteID) %>% 
  summarise(falltemp = mean3month[which(yday == 244)],
            wintertemp = mean3month[which(yday == 335)],
            allwintertemp = mean5month[which(yday == 305)],
            allgdd = max(accumdegday),
            firstfrostgdd = min(accumdegday[which(yday > 200 & hardfrost == TRUE)]),
            firstfrostdoy = min(yday[which(yday > 200 & hardfrost == TRUE)]))

# weather variation by site/year
gdd %>% 
  group_by(year) %>% 
  summarise(meantemp = mean((tmax..deg.c. + tmin..deg.c.)/2)) %>% 
  pull(meantemp) %>% 
  summary()
gdd %>% 
  group_by(SiteID) %>% 
  summarise(meantemp = mean((tmax..deg.c. + tmin..deg.c.)/2)) %>% 
  pull(meantemp) %>% 
  summary()
gdd %>% 
  group_by(SiteID, year) %>% 
  summarise(meantemp = mean((tmax..deg.c. + tmin..deg.c.)/2)) %>% 
  pull(meantemp) %>% 
  summary()

library(lmerTest)
# annual gdd trend
tidy(lmer(allgdd ~ year + (1|SiteID), data = weather_trends))
# later frost DOY, more GDD before first frost, almost equal to total trend in annual GDD
tidy(lmer(firstfrostdoy ~ year + (1|SiteID), data = weather_trends))
tidy(lmer(firstfrostgdd ~ year + (1|SiteID), data = weather_trends))
# Warmer 5 month winter temperature
tidy(lmer(allwintertemp~ year + (1|SiteID), data = weather_trends))
# Warmer 3 month winter temperature
tidy(lmer(wintertemp ~ year + (1|SiteID), data = weather_trends))
# Fall temperature Sept/Oct/Nov increasing
tidy(lmer(falltemp ~ year + (1|SiteID), data = weather_trends))


# no correlation between 5-month winter mean minimum temperature and first frost date (0.001)
# 3-month fall and 3-month winter temperatures are correlated (0.50)
# 3-month fall temperature and first frost DOY correlated (0.27)

GGally::ggpairs(weather_trends[, c("falltemp", "wintertemp", "allwintertemp", "firstfrostdoy")])


# Spatial and temporal variation in weather ----
# Compare change over time

ggplot(weather_trends, aes(x = year, y = firstfrostdoy, group = SiteID)) +
  # geom_point(shape = 95, size = 5, alpha = .33) +
  geom_point(
    position = position_jitter(width = .2, height = .1, seed = 0),
    size = 3, alpha = .1)+
  geom_abline(slope = 0.126, intercept = 58.6, color = "blue", linewidth = 1) +
  # scale_y_date(date_breaks = "1 month", date_labels = "%b") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  ylab("Day of first fall frost (-2°C)") +
  xlab("Year") +
  ggtitle("First fall frost", subtitle = "Trend and variation across sites and years")

ggplot(weather_trends, aes(x = year, y = firstfrostgdd, group = SiteID)) +
  geom_point(shape = 95, size = 5, alpha = .33) +
  # geom_point(
  #   position = position_jitter(width = .2, height = .1, seed = 0),
  #   size = 3, alpha = .1)+
  geom_abline(slope = 9.67, intercept = -16546, color = "blue", linewidth = 1) +
  # scale_y_date(date_breaks = "1 month", date_labels = "%b") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  ylab("Cumulative degree-days (5/30°C thresholds)\nby first fall frost (-2°C)") +
  xlab("Year") +
  ggtitle("Degree-days by first fall frost", subtitle = "Trend and variation across sites and years")

ggplot(weather_trends, aes(x = year, y = allwintertemp, group = SiteID)) +
  geom_point(shape = 95, size = 5, alpha = .33) +
  geom_abline(slope = 0.0287, intercept = -56.3, color = "blue", linewidth = 1) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  ylab("Winter mean minimum temperature (°C, Nov-Mar)") +
  xlab("Year") +
  ggtitle("Winter temperature", subtitle = "Trend and variation across sites and years")

ggplot(weather_trends, aes(x = year, y = allgdd, group = SiteID)) +
  geom_point(shape = 95, size = 5, alpha = .33) +
  geom_abline(slope = 9.85, intercept = -16840, color = "blue", linewidth = 1) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  ylab("Cumulative degree-days (5/30°C thresholds)") +
  xlab("Year") +
  ggtitle("Season length in degree-days", subtitle = "Trend and variation across sites and years")
  


# Imputation by generation ----

impdat <- genpop %>% 
  filter(CommonName %in% incl_species) %>% 
  filter(nyr >= 5, ObsSurvTotal >= 10, PosObsWeeks >= 3, trap.N.allgens >= 10, YearTotal >= 3) %>% 
  mutate(SiteYear = paste(SiteID, Year, sep = "_")) %>% 
  group_by(CommonName) %>% 
  mutate(maxbrood = base::max(as.numeric(gen))) %>% 
  group_by(CommonName, SiteYear, gen) %>%
  mutate(obs_vs_imp_weeks = weeks.obs.gen / (weeks.obs.gen + weeks.imp.gen),
         obs_vs_imp_total = total.obs.gen / (total.obs.gen + total.imp.gen))
ggplot(impdat %>% filter(gen == maxbrood, trap.N >=1), aes(x = obs_vs_imp_total))+
  geom_histogram() +
  facet_wrap(~CommonName, scales = "free")

ggplot(impdat %>% filter(gen == 1, trap.N >= 1), aes(x = obs_vs_imp_total))+
  geom_histogram() +
  facet_wrap(~CommonName, scales = "free")
