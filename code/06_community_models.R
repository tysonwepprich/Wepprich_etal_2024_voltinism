## Header ---------------------------
## Script name: 06_community_models.R
## Purpose of script: 
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## 
## ---

# Zuur model selection
# using moddat/test data from 07_last_generation_size.R
require(lmerTest)
library(ggpubr)

test <- moddat %>% 
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
    zlastann = lastloglambda - mean(lastloglambda)) %>% 
  group_by(SiteID) %>% 
  mutate(zfrostsite = mean(zfrost),
         zfrostann = zfrost - mean(zfrost),
         zwintersite = mean(zwinter, na.rm = TRUE),
         zwinterann = zwinter - mean(zwinter, na.rm = TRUE))

# 1. Beyond optimal model with all fixed effects
# 2. Fit increasing random effects complexity with REML
# 3. Compare nested fixed effects models with ML
# 4. Present final model with REML
modfixed <- glm(response ~ zdensann + (zyear + zordspec + zordsite + zordann)^3 + offset(off), family = poisson(link = "log"), data = test)


# is there a difference in estimates from glmer vs lmer? glmer is slower at fitting models but includes count size info
# results clearer with lmer interactions graphs

mod1 <- glmer(response ~ (zdensann + zyear + zordspec + zordsite + zordann)^2 + offset(off) + 
                # (1|SiteID) + (1|Year) +
                # (1 + (zyear + zdensann + zordsite + zordann)^2| CommonName) + 
                (1|CommonName) +
                (1|SiteID:CommonName) +
                (1|SpSiteYr),
              family = poisson(link = "log"), data = test, 
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

null_lam <- lmer(lastloglambda ~ zdensann + zyear + 
                   (1|CommonName) +
                   (1|SiteID:CommonName),
                 data = test, 
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

mod_lam1 <- lmer(lastloglambda ~ zdensann + zyear + (zordspec + zordsite + zordann)^3 + 
                   # (1|SiteID) +
                   # (1|Year) +
                   # (1 + zdensann + zyear + zordsite + zordann| CommonName) +
                   (1|CommonName) +
                   (1|SiteID:CommonName),
                 data = test, 
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

mod_lam2 <- lmer(lastloglambda ~ (zdensann + zyear + zordspec + zordsite + zordann)^2 + 
                   # (1|SiteID) +
                   # (1|Year) +
                   # (1 + zdensann + zyear + zordsite + zordann| CommonName) +
                   (1|CommonName) +
                   (1|SiteID:CommonName),
                 data = test, 
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

# photoperiod does worse than ordinal date! still!
# mod_lam2_photo <- lmer(lastloglambda ~ (zdensann + zyear + zphotospec + zphotosite + zphotoann)^3 + 
#                    # (1|SiteID) +
#                    # (1|Year) +
#                    # (1 + zdensann + zyear + zordsite + zordann| CommonName) +
#                    (1|CommonName) +
#                    (1|SiteID:CommonName),
#                  data = test, 
#                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

# this seems most in line with Gelman & Hill interaction with group-level predictor (zordspec)
mod_lam3 <- lmer(lastloglambda ~ zdensann + zyear + zordspec * (zordsite + zordann + zordsite:zordann) + 
                   (1|SiteID) +
                   # (1|Year) +
                   # (1 + zdensann + zyear + zordsite + zordann + zordsite:zordann| CommonName) +
                   (1|CommonName),
                 # (1|SiteID:CommonName),
                 data = test, 
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

saveRDS(mod_lam3, "LGmod_varyingslopes.rds")

mod_lam3a <- lmer(lastloglambda ~ zdensann + zyear + zordspec * (zordsite + zordann + zordsite:zordann) + 
                    # (1|SiteID) +
                    # (1|Year) +
                    # (1 + zdensann + zyear + zordsite + zordann + zordsite:zordann| CommonName) +
                    (1|CommonName) +
                    (1|SiteID:CommonName),
                  data = test, 
                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

saveRDS(mod_lam3a, "LGmod_varyingintcs.rds")
mod_lam3a <- readRDS("LGmod_varyingintcs.rds")

(step_res <- step(mod_lam3))
final <- get_model(step_res)
anova(final)
summary(final)

pltdat <- ggpredict(mod_lam3a, c("zordann", "zordsite[-.23,.23]", "zordspec[-1,1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Peak date annual variation", "Site mean phenology", "Species mean phenology")

# +/- 1 SD in species phenology is July 2, August 31
# 1 SD in scaled zord is approximately 30 days
# 1 SD in scaled zordsite is approximately 7 days (.23 * 30)
# 1 SD in scaled zordann is also approximately 7 days (.235 * 30)


plot(pltdat) + 
  labs(
    x = "Penultimate generation peak date (annual variation)",
    y = "Last generation size\n(population growth from penultimate generation)",
    title = NULL) +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Mean site\npeak date"))




tab_model(null_lam, mod_lam3a)


mod_lam5 <- lmer(lastloglambda ~ zyear + (zdensann + zordspec + zordsite + zordann)^3 + 
                   # (1|SiteID) +
                   # (1|Year) +
                   (1 + zyear + (zdensann + zordsite + zordann)^2| CommonName) +
                   # (1|CommonName) +
                   (1|SiteID:CommonName),
                 data = test, 
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

AIC(mod_lam3, mod_lam4, mod_lam5)


# OW lambda
test <- test[-which(is.na(test$zwinterann)),]
test <- test[-which(is.na(test$ow_lam)),]

mod2d <- lmer(ow_lam ~ (zordspec + zlastann + zdensann + zlastsite + zfrostann + zwinterann)^2 + 
                # (1 + zlastann * (zdensann + zlastsite + zfrostann + zwinterann) | CommonName) +
                (1|CommonName) + 
                (1|SiteID:CommonName), 
              data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

ggpredict(mod2d, c("zlastann", "zlastsite", "zfrostann", "zordspec")) %>% plot()


# quadratic effects for winter not exciting, but possible valley for zlastann that only adds .003 to the R2 when quadratic w/ interactions

# test <- test %>% filter(uniqSY >= 40 & uniqSYlam >= 100) # filter from 34 to 18 species

null_ow <- lmer(ow_lam ~ zdensann + zyear + 
                  (1|CommonName) +
                  (1|SiteID:CommonName), 
                data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# 
# 
# mod2d <- lmer(ow_lam ~ zdensann + zyear + (zordspec + poly(zlastann,2) + zlastsite + zfrostann + zwinterann)^2 + 
#                 # (1 + zlastann * (zdensann + zlastsite + zfrostann + zwinterann) | CommonName) +
#                 (1|CommonName) + 
#                 (1|SiteID:CommonName), 
#               data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# 
# ow <- lmer(ow_lam ~ zdensann + zyear + (zordspec + zlastann + zlastsite + zfrostann + zwinterann) + 
#                 # (1 + zlastann * (zdensann + zlastsite + zfrostann + zwinterann) | CommonName) +
#                 (1|CommonName) + 
#                 (1|SiteID:CommonName), 
#               data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# 
# ow2 <- lmer(ow_lam ~ zdensann + zyear + (zlastann + zlastsite + zfrostann + zwinterann)^2 + 
#              # (1 + zlastann * (zdensann + zlastsite + zfrostann + zwinterann) | CommonName) +
#              (1|CommonName) + 
#              (1|SiteID:CommonName), 
#            data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# 
# ow2a <- lmer(ow_lam ~ zdensann + zyear + (zordspec + zlastann + zlastsite + zfrostann + zwinterann)^2 +
#                (1 + zdensann + zyear + zlastann + zlastsite + zfrostann + zwinterann | CommonName) +
#                # (1|CommonName) + 
#                (1|SiteID:CommonName), 
#              data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

ow2b<- lmer(ow_lam ~ zdensann + zyear + zordspec * (zlastann + zlastsite + zfrostann + zwinterann)^2 + 
              # (1 + zdensann + zyear + (zlastann + zlastsite + zfrostann + zwinterann)^2 | CommonName) +
              (1|CommonName) +
              (1|SiteID:CommonName), 
            data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))


step_res <- step(ow2b, alpha.fixed = .1)
final <- get_model(step_res)
anova(final)
summary(final)


# SD zlastann = .92
# SD zlastsite = .99
# SD zordspec = .94
# SD zwinterann = .94
# SD zfrostann = .83



ggpredict(final, c("zlastann", "zlastsite", "zwinterann", "zordspec")) %>% plot()
ggpredict(final, c("zlastann", "zlastsite", "zfrostann", "zordspec")) %>% plot()
ggpredict(final, c("zlastann", "zwinterann", "zfrostann", "zordspec")) %>% plot()
ggpredict(final, terms = c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-.8,.8]", "zordspec[-1,1]")) %>% plot()


pltdat <- ggpredict(mod_lam3a, c("zordann", "zordsite[-.23,.23]", "zordspec[-1,1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Peak date annual variation", "Site mean phenology", "Species mean phenology")


tab_model(null_ow, final)

# 2 panels to save together Figure 3
pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-.8,.8]", "zordspec[-1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter temperature", "Species season")

pltdat1 <- plot(pltdat) + 
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Earlier season species", x = "", y = "") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Mean site\nlast generation\nsize"))

pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-.8,.8]", "zordspec[0]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter temperature", "Species season")

pltdat2 <- plot(pltdat) + 
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Mid season species",
       x = "",
       y = "Overwinter growth rate") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Mean site\nlast generation\nsize"))

pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-.8,.8]", "zordspec[1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter temperature", "Species season")

pltdat3 <- plot(pltdat) + 
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Later season species",
       x = "Last generation size (annual variation)",
       y = "") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Mean site\nlast generation\nsize"))


ggarrange(pltdat1, pltdat2, pltdat3, ncol = 1, labels = NULL, common.legend = TRUE, legend = "right")


# Frost x winter 2 panel
# 2 panels to save together Figure 3
pltdat <- ggpredict(final, c("zlastann[-3:3]", "zwinterann[-.95,.93]", "zfrostann[-.8,.8]", "zordspec[-1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Winter temperature", "Frost onset", "Species season")

pltdat1 <- plot(pltdat) + 
  scale_colour_brewer(name = "Winter\ntemperature", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Earlier season species", x = "", 
       y = "Overwinter growth rate") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Annual\nwinter\ntemperature"))

pltdat <- ggpredict(final, c("zlastann[-3:3]", "zwinterann[-.95,.93]", "zfrostann[-.8,.8]", "zordspec[1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Winter temperature", "Frost onset", "Species season")

pltdat2 <- plot(pltdat) + 
  scale_colour_brewer(name = "Winter\ntemperature", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Late season species",
       x = "Last generation size (annual variation)",
       y = "Overwinter growth rate") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Annual\nwinter\ntemperature"))

ggarrange(pltdat1, pltdat2, ncol = 1, labels = NULL, common.legend = TRUE, legend = "right")





# something like this might be best model, including 1 winter var at a time
# without varying slopes by species, will rely on species models to show variation
# seems like early season species have no impact of frost, but "trap" of larger LG only exists later in season and at colder sites.
ow3 <- lmer(ow_lam ~ zdensann + zyear + (zlastann + zlastsite + zfrostann + zwinterann)^3 + 
              # (1 + zlastann * (zdensann + zlastsite + zfrostann + zwinterann) | CommonName) +
              (1|CommonName) + 
              (1|SiteID:CommonName), 
            data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

ow3b <- lmer(ow_lam ~ zdensann + zyear + zordspec * (zlastann + zlastsite + zfrostann + zwinterann)^3 -zordspec:zlastann:zfrostann:zwinterann -
               zlastann:zfrostann:zwinterann -zordspec:zlastann:zlastsite:zfrostann - zlastann:zlastsite:zfrostann - zordspec:zlastsite:zfrostann:zwinterann -
               zordspec:zfrostann:zwinterann - zordspec:zlastsite:zfrostann -zordspec:zlastann:zlastsite:zwinterann - zordspec:zlastann:zlastsite -
               zordspec:zlastsite:zwinterann - zlastsite:zfrostann:zwinterann - zlastsite:zfrostann+ 
               # (1 + zlastann * (zdensann + zlastsite + zfrostann + zwinterann) | CommonName) +
               (1|CommonName) + 
               (1|SiteID:CommonName), 
             data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

# running into singularity issues when trying varying intercepts by species
# tried filtering rarer species, or reducing fixed effects (-zyear, one winter var at a time)
# not working out, but wanted to do it to match LG-lam models
#
ow2a<- lmer(ow_lam ~ zdensann + zyear + zordspec * (zlastann + zlastsite + zfrostann)^2 + 
              # (1 + zdensann + zyear + zlastann + zlastsite | CommonName) +
              (1|CommonName) +
              (1|SiteID), 
            data = test, REML = TRUE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))

ow2b <- lmer(ow_lam ~ zdensann + zordspec * (zlastann + zlastsite + zfrostann)^2 + 
               (1 + zdensann + (zlastann + zlastsite + zfrostann)^2   | CommonName) +
               # (1|CommonName) +
               (1|Year) +
               (1|SiteID), 
             data = test, REML = TRUE, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))



# trying to visualize differences between species
aaa <- test %>% ungroup %>% select(CommonName, zlastspec, zordspec) %>% distinct()
plot(aaa[,3:2])
text(aaa$zordspec, aaa$zlastspec, labels=aaa$CommonName)

