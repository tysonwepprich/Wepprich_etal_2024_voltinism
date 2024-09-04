## Header ---------------------------
## Script name: 06_community_models.R
## Purpose of script: Main analysis of community models (last generation size, overwinter growth rate)
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## 
## ---
source('code/01_data_prep.R')

library(lme4)
library(merTools)
library(lmerTest) # lmer and step functions masked
library(ggeffects)
library(broom)
library(broom.mixed)
library(ggpubr)
library(sjPlot)
library(performance)
moddat <- readRDS("data/modeling_data.rds")

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
    zlastann = lastloglambda - mean(lastloglambda)) %>% 
  group_by(SiteID) %>% 
  mutate(zfrostsite = mean(zfrost),
         zfrostann = zfrost - mean(zfrost),
         zwintersite = mean(zwinter, na.rm = TRUE),
         zwinterann = zwinter - mean(zwinter, na.rm = TRUE))


# Last generation size model ----
null_lam <- lmer(lastloglambda ~ zdensann + zyear + 
                   (1|CommonName) +
                   (1|SiteID:CommonName),
                 data = test, 
                 # weights = imputation_weights,
                 control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# check_model(null_lam)

mod_lam <- lmer(lastloglambda ~ zdensann + zyear + zordspec * (zordsite + zordann + zordsite:zordann) + 
                    (1|CommonName) +
                    (1|SiteID:CommonName),
                  data = test, 
                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# check_model(mod_lam)

mod_lam_wgt <- lmer(lastloglambda ~ zdensann + zyear + zordspec * (zordsite + zordann + zordsite:zordann) + 
                  (1|CommonName) +
                  (1|SiteID:CommonName),
                data = test, 
                weights = last_imputation_weights,
                control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# check_model(mod_lam_wgt)


(step_res <- step(mod_lam))
final <- get_model(step_res)
# anova(final)

(step_res <- step(mod_lam_wgt))
final_wgt <- get_model(step_res)

summary(final)
summary(final_wgt)

# Fig S2: Last generation response (community model) ----
pltdat <- ggpredict(final, c("zordann", "zordsite[-.23,.23]", "zordspec[-1,1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Peak date annual variation", "Site mean phenology", "Species mean phenology")

# +/- 1 SD in species phenology is July 2, August 31
# 1 SD in scaled zord is approximately 30 days
# 1 SD in scaled zordsite is approximately 7 days (.23 * 30)
# 1 SD in scaled zordann is also approximately 7 days (.235 * 30)


pltdat <- plot(pltdat) + 
  labs(
    x = "Annual variation in phenology (standard deviation)",
    y = "Last generation size\nlog(last generation / penultimate generation)",
    title = NULL) +
  scale_colour_brewer(palette = "Set1", labels = c("-1 SD (earlier, warmer site)", "+1 SD (later, cooler site)")) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw(base_size = 18) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(legend.position = c(.15, .85)) +
  guides(color=guide_legend(title="Site mean\nphenology"))
pltdat
ggsave(filename = "figS2.tif", path = "figures", device='tiff', dpi=600)


# Table S1: Last generation size ----
tab_model(null_lam, final, show.se = TRUE, show.ci = FALSE, show.icc = FALSE, 
          show.stat = TRUE, digits = 3, digits.re = 3, CSS = list(css.tdata = '+padding:0.1cm;')
)
# Table S6: Imputation robustness check ----
tab_model(final, final_wgt, show.se = TRUE, show.ci = FALSE, show.icc = FALSE, 
          show.stat = TRUE, digits = 3, digits.re = 3, CSS = list(css.tdata = '+padding:0.1cm;'))


# Overwinter growth rate model ----
test <- test[-which(is.na(test$zwinterann)),]
test <- test[-which(is.na(test$ow_lam)),]

null_ow <- lmer(ow_lam ~ zdensann + zyear + 
                  (1|CommonName) +
                  (1|SiteID:CommonName), 
                data = test, control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
# check_model(null_ow)

mod_ow <- lmer(ow_lam ~ zdensann + zyear + zordspec * (zlastann + zlastsite + zfrostann + zwinterann)^2 + 
                     (1|CommonName) +
                     (1|SiteID:CommonName), 
                   data = test, 
                   control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
# check_model(mod_ow)

mod_ow_wgt <- lmer(ow_lam ~ zdensann + zyear + zordspec * (zlastann + zlastsite + zfrostann + zwinterann)^2 + 
              (1|CommonName) +
              (1|SiteID:CommonName), 
            data = test, 
            weights = last_imputation_weights + first_imputation_weights,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e7)))
# check_model(mod_ow_wgt)

step_res <- step(mod_ow, alpha.fixed = .1)
final <- get_model(step_res)
# anova(final)


step_res <- step(mod_ow_wgt, alpha.fixed = .1)
final_wgt <- get_model(step_res)

summary(final)
summary(final_wgt)

# Table S2: Overwinter population growth rate ----
tab_model(null_ow, final,  show.se = TRUE, show.ci = FALSE, show.icc = FALSE, 
          show.stat = TRUE, digits = 3, digits.re = 3, CSS = list(css.tdata = '+padding:0.1cm;'))

# Table S7: Imputation robustness check ----
tab_model(final, final_wgt,  show.se = TRUE, show.ci = FALSE, show.icc = FALSE, 
          show.stat = TRUE, digits = 3, digits.re = 3, CSS = list(css.tdata = '+padding:0.1cm;'))




# ggpredict(final, c("zlastann", "zlastsite", "zwinterann", "zordspec")) %>% plot()
ggpredict(final, c("zlastann", "zlastsite", "zfrostann", "zordspec")) %>% plot()
# ggpredict(final, c("zlastann", "zwinterann", "zfrostann", "zordspec")) %>% plot()
# ggpredict(final, terms = c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-.8,.8]", "zordspec[-1,1]")) %>% plot()



# Fig 4: Overwinter: winter onset ----
# 2 panels to save together Figure 3
pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zfrostann[-1,1]", "zordspec[-1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter onset", "Species season")

pltdat1 <- plot(pltdat) + 
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  scale_y_continuous(limits = c(-1.9, .6), breaks = c(-1.5, -1.0, -0.5, 0, 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = NULL, x = NULL, y = "") +
  annotate("text", x = -1.72, y = .38, label = "Earlier species\nJuly 2 peak", size = 5) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(legend.position = c(.15, .2)) +
  guides(color=guide_legend(title="Site mean\nlast generation size"))
# pltdat1

pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zfrostann[-1,1]", "zordspec[0]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter onset", "Species season")

pltdat2 <- plot(pltdat) + 
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  scale_y_continuous(limits = c(-1.9, .6), breaks = c(-1.5, -1.0, -0.5, 0, 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = NULL,
       x = NULL,
       y = "") +
  annotate("text", x = -1.72, y = .38, label = "Mid species\nAugust 1 peak", size = 5) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(legend.position = "none") +
  guides(color=guide_legend(title="Site mean\nLG size"))
# pltdat2

pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zfrostann[-1,1]", "zordspec[1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter onset", "Species season")

pltdat3 <- plot(pltdat) + 
  scale_y_continuous(limits = c(-1.9, .6), breaks = c(-1.5, -1.0, -0.5, 0, 0.5)) +
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = NULL,
       x = NULL,
       y = "") +
  annotate("text", x = -1.72, y = .38, label = "Later species\nAugust 31 peak", size = 5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Site mean\nLG size"))
# pltdat3

ggarrange(pltdat1, pltdat2, pltdat3, ncol = 1, labels = NULL, common.legend = FALSE)
ggsave(filename = "fig4.tif", path = "figures", device='tiff', dpi=600, width = 9, height = 12.5, units = "in")
#900x1200



# Fig S4: Overwinter: winter temperature ----
# 2 panels to save together Figure S4
pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-1,1]", "zordspec[-1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter temperature", "Species season")

pltdat1 <- plot(pltdat) + 
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  scale_y_continuous(limits = c(-2, .6), breaks = c(-1.5, -1.0, -0.5, 0, 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = NULL, x = NULL, y = "") +
  annotate("text", x = -1.72, y = .38, label = "Earlier species\nJuly 2 peak", size = 5) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(legend.position = c(.15, .2)) +
  guides(color=guide_legend(title="Site mean\nlast generation size"))
# pltdat1

pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-1,1]", "zordspec[0]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter temperature", "Species season")

pltdat2 <- plot(pltdat) + 
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  scale_y_continuous(limits = c(-2, .6), breaks = c(-1.5, -1.0, -0.5, 0, 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(title = NULL,
       x = NULL,
       y = "") +
  annotate("text", x = -1.72, y = .38, label = "Mid species\nAugust 1 peak", size = 5) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  theme(legend.position = "none") +
  guides(color=guide_legend(title="Site mean\nLG size"))
# pltdat2

pltdat <- ggpredict(final, c("zlastann[-3:3]", "zlastsite[-1,1]", "zwinterann[-1,1]", "zordspec[1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Site mean last generation size", "Winter temperature", "Species season")

pltdat3 <- plot(pltdat) + 
  scale_y_continuous(limits = c(-2, .6), breaks = c(-1.5, -1.0, -0.5, 0, 0.5)) +
  scale_colour_brewer(name = "Site mean\nLG size", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = NULL,
       x = NULL,
       y = "") +
  annotate("text", x = -1.72, y = .38, label = "Later species\nAugust 31 peak", size = 5) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Site mean\nLG size"))
# pltdat3

ggarrange(pltdat1, pltdat2, pltdat3, ncol = 1, labels = NULL, common.legend = FALSE)
ggsave(filename = "figS4.tif", path = "figures", device='tiff', dpi=600,  width = 9, height = 12.5, units = "in")


# Fig S5: Frost x temperature -----
pltdat <- ggpredict(final, c("zlastann[-3:3]", "zwinterann[-1,1]", "zfrostann[-1,1]", "zordspec[-1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Winter temperature", "Winter onset", "Species season")

pltdat1 <- plot(pltdat) + 
  scale_colour_brewer(name = "Winter\ntemperature", palette = "Set1", direction = -1, labels = c("-1 SD (cooler)", "+1 SD (warmer)")) +
  scale_fill_brewer(palette = "Set1", direction = -1) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = "Earlier season species", x = "", 
       y = "Overwinter growth rate") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  guides(color=guide_legend(title="Annual\nwinter\ntemperature"))

pltdat <- ggpredict(final, c("zlastann[-3:3]", "zwinterann[-1,1]", "zfrostann[-1,1]", "zordspec[1]")) 
# change attributes so facets labeled (3rd one only?)
attr(pltdat, "terms") <- c("Last generation size", "Winter temperature", "Winter onset", "Species season")

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
ggsave(filename = "figS5.tif", path = "figures", device='tiff', dpi=600,  width = 9, height = 9, units = "in")

