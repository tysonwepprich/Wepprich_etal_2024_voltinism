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


genpop <- readRDS("data/genpops.rds")

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
         lastloglambda = log((trap.N[maxbrood]+1) / (trap.N[maxbrood-1]+1)))


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

# 
# # Cool plot of mismatch by GDD left versus first hard frost. Generally, year is more variable than site, 
# # although in some years there is a bifurcation among sites for first frost timing. 
# # Imagine some years would be good for all species attempting a last generation.
# ggplot(gdd_left %>% filter(yday == 225), aes(x = gddmismatch, y = gddmismatchfrost, color = as.factor(year))) +
#   geom_point() + facet_wrap(~year) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal()
# ggplot(gdd_left %>% group_by(year, yday) %>% summarise(meangdd = mean(accumdegday)), aes(x = yday, y = meangdd, group = year, color = year)) + 
#   geom_line()
# ggplot(gdd_left %>% filter(yday == 365), aes(x = year, y = accumdegday)) + 
#   geom_point() + geom_smooth(method = "lm") + xlab("Year") + ylab("Annual growing degree-days (base 5°C)")

# Covariates at brood peaks ----
# UNCOMMENT FOR FIRST RUN

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
#   # mean temperature/precip near peak, larger temperature window makes it too collinear with phenology changes (ordinal date)
#   # don't use these after all
#   tempwindow <- 15
#   meantemp <- temp %>%
#     filter(yday >= ord - tempwindow & yday <= ord + tempwindow)
#   meantemp <- mean(meantemp$daymeantemp)
#   meanprec <- temp %>%
#     filter(yday >= ord - tempwindow & yday <= ord + tempwindow)
#   meanprec <- mean(meanprec$prcp..mm.day.)
#   tempwindow <- 7.5
#   meantemphalf <- temp %>%
#     filter(yday >= ord - tempwindow & yday <= ord + tempwindow)
#   meantemphalf <- mean(meantemphalf$daymeantemp)
#   meanprechalf <- temp %>%
#     filter(yday >= ord - tempwindow & yday <= ord + tempwindow)
#   meanprechalf <- mean(meanprechalf$prcp..mm.day.)
# 
#   return(data.frame(ord,
#              photo,
#              gddleft,
#              gddexpected,
#              gddmismatch,
#              gddleftfrost,
#              gddexpectedfrost,
#              gddmismatchfrost,
#              meantemp,
#              meanprec,
#              meantemphalf,
#              meanprechalf))
# }
# 
# 
# system.time({covs <- popdat %>%
#   group_by(CommonName, SiteYear, gen) %>%
#   do(gdd2vars(.$SiteYear, .$trap.mu.gdd))
# })
# saveRDS(covs, "covs.rds")

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
         firstYear = lead(YearNum, 1)) %>% 
  filter(firstYear == (YearNum + 1)) %>% 
  filter(firstN > 1) %>% 
  dplyr::select(SiteID, Year, SiteYear, CommonName, firstN:firstYear)

moddat <- moddat %>% 
  left_join(firstdat)

# Winter covariates ----
# Used 5 month winter mean minimum, but labeled 6 month accidentally
winter <- gdd %>% 
  group_by(SiteID) %>% 
  arrange(SiteDate) %>% 
  mutate(mean3month = zoo::rollmean(tmin..deg.c., 90, fill = NA, align = "left"),
         mean6month = zoo::rollmean(tmin..deg.c., 150, fill = NA, align = "left")) %>% 
  group_by(SiteYear, year, SiteID) %>% 
  summarise(falltemp = mean3month[which(yday == 244)],
            wintertemp = mean3month[which(yday == 365)],
            allwintertemp = mean6month[which(yday == 305)],
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

# # plot how some years have many species producing extra generation (2010-2012 for example)
# moddat <- readRDS("data/modeling_data.rds")
# pltpop <- moddat %>% 
#   group_by(CommonName, SiteID) %>% 
#   mutate(sitemean = mean(lastprop),
#          siteyearmean = lastprop - sitemean)
# ggplot(pltpop %>% filter(CommonName == "Cabbage White") , aes(x = sitemean, y = siteyearmean, color = Year)) +
#   geom_point() + facet_wrap(~Year) # + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal()
# ggplot(pltpop %>% filter(CommonName == "Pearl Crescent"), aes(x = YearNum, y = siteyearmean)) +
#   geom_point() + geom_smooth() + facet_wrap(~CommonName) # + geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal()
# 
# 

# Fig 2D ----
# example species
lgvar <- moddat %>%
  filter(CommonName == "Pearl Crescent", year %in% c(2008:2017), YearTotal >= 50, nyr >= 15, photo > 12.2) %>%
  mutate(yrgroup = ifelse(year %in% c(2010, 2011, 2012, 2016, 2017), "Warm years", "Cool years"),
         sitegroup = ifelse(region %in% c("NE", "NW"), "Cool sites", "Warm sites"))
lgvarplt <- ggplot(data = lgvar, aes(x = photo,
                                      y = lastloglambda,
                                      color = yrgroup))+
  geom_point()+
  scale_x_reverse() +
  scale_color_brewer(name = NULL, palette = "Set1", direction = -1) +
  scale_fill_brewer(name = NULL, palette = "Set1", direction = -1) +
  facet_wrap(~sitegroup, ncol = 1) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = c(.25, .1)) +
  labs(
    x = "Day length at penultimate peak",
    y = "Last generation size",
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
         mean6month = zoo::rollmean((tmin..deg.c. + tmax..deg.c.)/2, 180, fill = NA, align = "left")) %>% 
  group_by(SiteYear, year, SiteID) %>% 
  summarise(falltemp = mean3month[which(yday == 244)],
            wintertemp = mean3month[which(yday == 335)],
            allwintertemp = mean6month[which(yday == 305)],
            siteyrtotalgdd = max(accumdegday),
            firstfrostgdd = min(accumdegday[which(yday > 200 & hardfrost == TRUE)]),
            firstfrostdoy = min(yday[which(yday > 200 & hardfrost == TRUE)]),
            allgdd = accumdegday[which(yday == 365)])

# annual gdd trend
tidy(lmer(allgdd ~ year + (1|SiteID), data = weather_trends))
# later frost DOY, more GDD before first frost, almost equal to total trend in annual GDD
tidy(lmer(firstfrostdoy ~ year + (1|SiteID), data = weather_trends))
tidy(lmer(firstfrostgdd ~ year + (1|SiteID), data = weather_trends))
# Warmer 6 month winter temperature
tidy(lmer(allwintertemp~ year + (1|SiteID), data = weather_trends))
# Warmer 3 month winter temperature
tidy(lmer(wintertemp ~ year + (1|SiteID), data = weather_trends))
# Fall temperature Sept/Oct/Nov increasing
tidy(lmer(falltemp ~ year + (1|SiteID), data = weather_trends))


