## Header ---------------------------
## Script name: 07_voltinism_traits.R
## Purpose of script: Analysis of species' traits correlated with population trends, phylogenetic analysis
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## 
## ---
library(corrplot)
library(ggrepel)
library(ape)
library(caper)
library(lme4)
library(lmerTest)
library(broom)
library(ggrepel)

source('code/01_data_prep.R')
moddat <- readRDS("data/modeling_data.rds")

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


# Methods verification
# shows how GDD always had lower class uncertainty that DOY
mmcompare <- readRDS("data/mixmodcomparison.rds")
mmused <- mmcompare %>% 
  filter(uncertainty > 0.000000001) %>% 
  # filter(time != "doy") %>%
  group_by(species) %>% 
  # mutate(maxgen = max(gen)) %>% 
  # filter(gen == maxgen) %>%
  filter(uncertainty == min(uncertainty))

species_table <- mmused %>% 
  dplyr::select(species, gen, prop) %>% 
  group_by(species) %>% 
  summarise(Generations = paste(paste0(as.character(c(round(prop*100))), "%"), collapse = "/")) %>% 
  ungroup() %>% 
  mutate(CommonName = species) %>% 
  dplyr::select(-species)

species_table <- left_join(species_table, moddat %>% select(CommonName, uniqSY, uniqSYlam) %>% distinct() %>%  arrange(uniqSY) %>% data.frame())

penult_doy <- moddat %>% 
  group_by(CommonName) %>% 
  summarise(meandoy = round(median(ord, na.rm = TRUE))) %>% 
  mutate(meandoy = format(strptime(paste(meandoy, "2000", sep = "-"), "%j-%Y"), "%d-%b"))          

species_table <- left_join(species_table, penult_doy)


# Look at species parameters from 05_species_models.R

allpars <- readRDS("data/species_models.rds")

# simulated benefit of last generation
# 8 species not enough data unless filter changed from 10/site to 5/site
# 8 negative, only EuroSkip & Lil Glassywing signif.
# 22 positive, 14 signif
lg_sim1 <- allpars %>% filter(model == "lg_sim1")
lg_sim1 %>% filter(!is.na(std.error)) %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)

# local adaptation
allpars %>% filter(model == "vdp_2") %>% 
  filter(term == "zordsite") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)
# pheno plasticity
allpars %>% filter(model == "vdp_2") %>% 
  filter(term == "zordann") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)
# signif co-gradient: 7 species
sp1 <- allpars %>% filter(model == "vdp_2") %>% 
  filter(term == "zordann") %>% filter(p.value < 0.05) %>%  arrange(estimate) %>%  data.frame() %>% round_df(digits = 2) %>%  pull(CommonName)
sp2 <- allpars %>% filter(model == "vdp_2") %>% 
  filter(term == "zordsite") %>% filter(p.value < 0.05) %>%  arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)%>%  pull(CommonName)
allpars %>% filter(model == "vdp_2", CommonName %in% sp1[sp1 %in% sp2]) %>% data.frame() %>%  round_df(digits = 2)



# ow model
allpars %>% filter(model == "ow") %>% 
  filter(term == "zdensann") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)

# trend
# 20 negative/10 positive (or 6 sig neg/3 sig pos)
allpars %>% filter(model == "trend") %>% 
  filter(term == "zyear") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 3)



# compare some traits with trends
traits <- allpars %>% 
  group_by(CommonName) %>% 
  summarise(LG_trend = estimate[which(term == "zyear" & model == "lastgen")],
            local_adapt = estimate[which(term == "zordsite" & model == "vdp_2")],
            annual_plast = estimate[which(term == "zordann" & model == "vdp_2")],
            trend = estimate[which(term == "zyear" & model == "trend")],
            lg_sim = ifelse(is.na(std.error[which(term == "(Intercept)" & model == "lg_sim1")]), NA, estimate[which(term == "(Intercept)" & model == "lg_sim1")]))
traits <- left_join(traits, distinct(moddat2[,c("CommonName", "maxbrood", "zordspec", "zlastspec")]))


# Species trait table ----
species <- read.csv("data/species_names.csv") %>% 
  dplyr::select(CommonName, Genus, Species) %>%
  distinct() %>% 
  mutate(Latin = paste(Genus, Species, sep = " ")) 
species_table <- left_join(species_table, species)

# # LG Trend
# lgtrend_trait <- bind_rows(outlist) %>% filter(term == "zyear") %>% 
#   select(estimate, std.error, CommonName) %>% 
#   mutate(LG_trend = paste0(round(estimate,3), "  [", round(estimate-1.96*std.error,3), ", ", round(estimate+1.96*std.error,3), "]"))
# species_table <- left_join(species_table, lgtrend_trait[,c("CommonName", "LG_trend")])

modtraits <- allpars %>% 
  group_by(CommonName) %>% 
  summarise(LG_trend = estimate[which(term == "zyear" & model == "lastgen")],
            LG_trend_se = std.error[which(term == "zyear" & model == "lastgen")],
            local_adapt = estimate[which(term == "zordsite" & model == "vdp_2")],
            local_adapt_se = std.error[which(term == "zordsite" & model == "vdp_2")],
            annual_plast = estimate[which(term == "zordann" & model == "vdp_2")],
            annual_plast_se = std.error[which(term == "zordann" & model == "vdp_2")],
            trend = estimate[which(term == "zyear" & model == "trend")],
            trend_se = std.error[which(term == "zyear" & model == "trend")],
            lg_sim = estimate[which(term == "(Intercept)" & model == "lg_sim1")],
            lg_sim_se = std.error[which(term == "(Intercept)" & model == "lg_sim1")]) %>% 
  mutate(LG_trend = paste0(round(LG_trend,3), "  [", round(LG_trend-1.96*LG_trend_se,3), ", ", round(LG_trend+1.96*LG_trend_se,3), "]"),
         LGspatial = paste0(round(local_adapt,3), "  [", round(local_adapt-1.96*local_adapt_se,3), ", ", round(local_adapt+1.96*local_adapt_se,3), "]"),
         LGannual = paste0(round(annual_plast,3), "  [", round(annual_plast-1.96*annual_plast_se,3), ", ", round(annual_plast+1.96*annual_plast_se,3), "]"),
         LGsim = paste0(round(lg_sim,3), "  [", round(lg_sim-1.96*lg_sim_se,3), ", ", round(lg_sim+1.96*lg_sim_se,3), "]"),
         PopTrend = paste0(round(trend,3), "  [", round(trend-1.96*trend_se,3), ", ", round(trend+1.96*trend_se,3), "]"))
species_table <- left_join(species_table, modtraits[,c("CommonName", "LG_trend", "LGspatial", "LGannual", "LGsim", "PopTrend")])
species_table$Sample <- paste(species_table$uniqSY, species_table$uniqSYlam, sep = "/")

species_table_out <- species_table[which(!is.na(species_table$LG_trend)), c("Latin", "Sample", "Generations", "meandoy", "LG_trend", "LGspatial", "LGannual", "LGsim", "PopTrend")] %>% 
  arrange(Latin)

# not included in manuscript
write.csv(species_table_out, "Table1.csv", row.names = FALSE)

# Fig. 4: caterpillar plot ----

modtraits <- allpars %>% 
  group_by(CommonName) %>% 
  summarise(LG_trend = estimate[which(term == "zyear" & model == "lastgen")],
            LG_trend_se = std.error[which(term == "zyear" & model == "lastgen")],
            local_adapt = estimate[which(term == "zordsite" & model == "vdp_2")],
            local_adapt_se = std.error[which(term == "zordsite" & model == "vdp_2")],
            annual_plast = estimate[which(term == "zordann" & model == "vdp_2")],
            annual_plast_se = std.error[which(term == "zordann" & model == "vdp_2")],
            trend = estimate[which(term == "zyear" & model == "trend")],
            trend_se = std.error[which(term == "zyear" & model == "trend")],
            lg_sim = estimate[which(term == "(Intercept)" & model == "lg_sim1")],
            lg_sim_se = std.error[which(term == "(Intercept)" & model == "lg_sim1")])

catdat1 <- modtraits %>% 
  dplyr::select(-ends_with("_se")) %>% 
  pivot_longer(cols = LG_trend:lg_sim, names_to = "Trait", values_to = "Estimate")
catdat2 <- modtraits %>% 
  dplyr::select(ends_with("_se")) %>% 
  pivot_longer(cols = LG_trend_se:lg_sim_se, names_to = "Trait_SE", values_to = "SE")
catdat <- bind_cols(catdat1, catdat2) %>% 
  mutate(upr = Estimate + 1.96*SE,
         lwr = Estimate - 1.96*SE) %>%
  mutate(pnt_shape = ifelse(upr*lwr > 0, "sig", "nonsig")) %>% 
  left_join(species)

sp_levels <- catdat %>% filter(Trait == "trend") %>% arrange(-Estimate) %>% pull(Latin) %>% unique()
catdat$Trait <- factor(catdat$Trait, levels = c("LG_trend", "local_adapt", "annual_plast", "lg_sim", "trend"), ordered = TRUE)
levels(catdat$Trait) <- c("Trend in\nLG size", "LG response to\nsite phenology", "LG response to\nannual phenology", 
                                  "Simulated effect\nof larger LG","Population\ntrend")
catdat$Latin_ordered = factor(catdat$Latin, levels=sp_levels, ordered = TRUE)

ggplot(catdat, aes(y = Latin_ordered, x = Estimate, color = pnt_shape)) + 
  geom_point() + 
  scale_color_manual(name = NULL, values = c("gray", "black")) +
  geom_pointrange(aes(xmin = lwr, 
                      xmax = upr)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = seq(5.5, 25.5,length.out = 5), alpha = .1) +
  xlab("Estimate (95% CI)") +
  ylab(NULL) +
  facet_grid(~Trait, scales = "free") +
  theme(legend.position = "none") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(face = "italic"))
ggsave(filename = "fig4.tif", path = "figures", device='tiff', dpi=600)


species <- read.csv("data/species_names.csv") %>% 
  dplyr::select(CommonName, Genus, Species) %>%
  distinct() %>% 
  mutate(Latin = paste(Genus, Species, sep = " "),
         Latin_abbrev = paste0(str_sub(Genus, 1, 1), ". ", Species)) 
plttraits <- left_join(traits, species)

# Fig. 5B ----
ggplot(plttraits, aes(x = LG_trend, y = lg_sim, label = Latin_abbrev, color = trend)) +
  geom_point() +
  scale_color_viridis(name = "Population\ntrend", option = "C", begin = .8, end = 0) +
  geom_text_repel(size = 5, fontface = "italic") +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .4) +
  xlab("Observed trend in last generation size") +
  ylab("Simulated change in overwinter population growth\nwith 1 std. dev. larger last generation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = c(.875, .2))

ggsave(filename = "fig5B.tif", path = "figures", device='tiff', dpi=600)


# correlations without bootstrap ----
# not used in manuscript

traits_clean <- traits
names(traits_clean) <- c("CommonName", "Trend in\nLG size", "LG response to\nsite phenology", "LG response to\nannual phenology", "Population\ntrend",
                         "Simulated effect\nof larger LG", "Max # of\ngenerations", "Mean penult.\npeak date",
                         "Mean\nLG size")
traits_clean <- traits_clean[,c(1, 8, 9, 2, 3,4,6,5)]
GGally::ggpairs(traits_clean[, -1])

M <- cor(as.matrix(traits_clean[,-1]))
testRes = cor.mtest(M, conf.level = 0.95)
corrplot(M, type = "upper", addCoef.col = 'black')
 



# Fig 5A correlations with bootstrap----
set.seed(1)
sp <- sort(unique(moddat2$CommonName))

cordat <- array(data = NA, dim = c(30, 7, 1000))
for (i in 1:30){
  
  tmpars <- allpars %>% filter(CommonName == sp[i])
  
  # species phenology
  cordat[i,1,] <- sur::boot.mean(moddat2$zord[moddat2$CommonName == sp[i]], 1000)$bootstrap.samples
  # species LG size
  cordat[i,2,] <- sur::boot.mean(moddat2$lastloglambda[moddat2$CommonName == sp[i]], 1000)$bootstrap.samples
  # species LG trend
  cordat[i,3,] <- rnorm(n = 1000, mean = tmpars$estimate[which(tmpars$term == "zyear" & tmpars$model == "lastgen")], sd = tmpars$std.error[which(tmpars$term == "zyear" & tmpars$model == "lastgen")])
  # spatial var
  cordat[i,4,] <- rnorm(n = 1000, mean = tmpars$estimate[which(tmpars$term == "zordsite" & tmpars$model == "vdp_2")], sd = tmpars$std.error[which(tmpars$term == "zordsite" & tmpars$model == "vdp_2")])
  # temporal var
  cordat[i,5,] <- rnorm(n = 1000, mean = tmpars$estimate[which(tmpars$term == "zordann" & tmpars$model == "vdp_2")], sd = tmpars$std.error[which(tmpars$term == "zordann" & tmpars$model == "vdp_2")])
  # simulation effect
  cordat[i,6,] <- rnorm(n = 1000, mean = tmpars$estimate[which(tmpars$term == "(Intercept)" & tmpars$model == "lg_sim1")], sd = tmpars$std.error[which(tmpars$term == "(Intercept)" & tmpars$model == "lg_sim1")])
  # statewide trend
  cordat[i,7,] <- rnorm(n = 1000, mean = tmpars$estimate[which(tmpars$term == "zyear" & tmpars$model == "trend")], sd = tmpars$std.error[which(tmpars$term == "zyear" & tmpars$model == "trend")])
}

M.bs <- apply(cordat, 3, FUN = cor) 
M.bs <- array(data = M.bs, dim = c(7,7,1000))
M.bs <- apply(M.bs, c(1,2), FUN = stats::quantile, probs = c(0.025, .5, 0.975))

Mcoef <- M.bs[2,,]
attr(Mcoef, "dimnames") <- attr(M, "dimnames")
Mlo <- M.bs[1,,]
attr(Mlo, "dimnames") <- attr(M, "dimnames")
Mhi <- M.bs[3,,]
attr(Mhi, "dimnames") <- attr(M, "dimnames")

tiff('figures/fig5A.tiff', units="in", width=7, height=7, res=600, compression = 'lzw')

corrplot::corrplot(Mcoef, type = "upper", addCoef.col = 'black', tl.col	= "black", col = COL2('BrBG', 100), number.digits = 2, diag = TRUE)

conf <- paste0("(", format(Mlo, digits=1, nsmall = 2), ",", format(Mhi, digits=1, nsmall = 2), ")")
xs <- row(Mlo)
ys <- (ncol(Mlo)+1) - col(Mlo)
text(xs[lower.tri(xs)], ys[lower.tri(xs)], conf[lower.tri(xs)], pos=1, cex=.75)
# bug in corrplot, can't add coef for bottom right corner, but only when diag = FALSE in corrplot function above!
paste0("(", format(Mlo[6,7], digits=1, nsmall = 2), ",", format(Mhi[6,7], digits=1, nsmall = 2), ")")

dev.off()


# visualization of trait plots with correlation bootstrap error
trait.bs <- apply(cordat, c(1,2), FUN = stats::quantile, probs = c(0.025, .5, 0.975))
trait.bs.lo <- data.frame(CommonName = sp, trait.bs[1,,], quant = "lo")
trait.bs.med <- data.frame(CommonName = sp, trait.bs[2,,], quant = "med")
trait.bs.hi <- data.frame(CommonName = sp, trait.bs[3,,], quant = "hi")
trait.df <- rbind(trait.bs.lo, trait.bs.med, trait.bs.hi)
names(trait.df) <- c(names(traits)[c(1,8,9,2,3,4,6,5)], "quant")

trend_wide <- trait.df %>% 
  dplyr::select(CommonName, trend, quant) %>% 
  pivot_wider(id_cols = CommonName, names_from = quant, values_from = trend)
names(trend_wide) <- c("CommonName", "trend_lo", "trend", "trend_hi")
lg_trend_wide <- trait.df %>% 
  dplyr::select(CommonName, LG_trend, quant) %>% 
  pivot_wider(id_cols = CommonName, names_from = quant, values_from = LG_trend)
names(lg_trend_wide) <- c("CommonName", "LG_trend_lo", "LG_trend", "LG_trend_hi")

plt_trait <- merge(trend_wide, lg_trend_wide)

ggplot(plt_trait, aes(x = LG_trend, y = trend, label = CommonName)) +
  geom_point() +
  geom_text_repel(size = 4) +
  geom_linerange(aes(ymin = trend_lo, ymax = trend_hi), alpha = .3) +
  geom_linerange(aes(xmin = LG_trend_lo, xmax = LG_trend_hi), alpha = .3) +
  # scale_y_continuous(limits=c(0,100), expand = expansion(mult = c(0, 0))) +
  # geom_hline(yintercept = meantrend, linetype = "dashed", alpha = .4) +
  # geom_vline(xintercept = meansim, linetype = "dashed", alpha = .4) +
  xlab("Observed trend in last generation size") +
  ylab("Observed statewide population trend (1st generation, 1996-2022)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

OW_sim_wide <- trait.df %>% 
  dplyr::select(CommonName, lg_sim, quant) %>% 
  pivot_wider(id_cols = CommonName, names_from = quant, values_from = lg_sim)
names(OW_sim_wide) <- c("CommonName", "LG_sim_lo", "LG_sim", "LG_sim_hi")
plt_trait <- merge(trend_wide, OW_sim_wide)

ggplot(plt_trait, aes(x = LG_sim, y = trend, label = CommonName)) +
  geom_point() +
  geom_text_repel(size = 4) +
  geom_linerange(aes(ymin = trend_lo, ymax = trend_hi), alpha = .3) +
  geom_linerange(aes(xmin = LG_sim_lo, xmax = LG_sim_hi), alpha = .3) +
  # scale_y_continuous(limits=c(0,100), expand = expansion(mult = c(0, 0))) +
  # geom_hline(yintercept = meantrend, linetype = "dashed", alpha = .4) +
  # geom_vline(xintercept = meansim, linetype = "dashed", alpha = .4) +
  xlab("Simulated change in overwinter population growth\nwith 1 std. dev. larger last generation") +
  ylab("Observed statewide population trend (1st generation, 1996-2022)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# phylogenetic without bootstrap ----
trends <- traits

best.tree <- read.nexus(file = "data/RAxML_bipartitions.OHbflyBOLD_JRA_cds_aligned.gBlocks.consensus.famsConstrained.nex")
tr <- root(best.tree, outgroup = 49, resolve.root = TRUE) #Bombyxmori outgroup, silk moth
tr$tip.label <- gsub(pattern = "'", replacement = "", x = tr$tip.label, fixed = TRUE)


species <- read.csv("data/species_names.csv") %>% 
  dplyr::select(CommonName, Genus, Species) %>%
  distinct() %>% 
  mutate(Latin = paste(Genus, Species, sep = " ")) 
trends <- left_join(trends, species)

spp <- read.csv("data/OHbflyID.csv")
trends <- trends %>% 
  left_join(distinct(spp[,c("CommonName", "species_name")])) %>% 
  mutate(species_name = gsub(pattern = " ", replacement = "", x = species_name, fixed = TRUE))
trends$species_name[which(trends$CommonName %in% c("Azures", "Spring/Summer Azure"))] <- "Celastrinaladon"

trends <- trends[which(trends$species_name %in% tr$tip.label), ]
trends <- droplevels(trends)
newtips <- trends$Latin[match(tr$tip.label, trends$species_name)]
tr$tip.label <- newtips


trend <- trends %>% pull(trend)
names(trend) <- trends %>% pull(Latin)
tr <- keep.tip(tr, tip = names(trend))

trends <- data.frame(trends)
row.names(trends) <- trends$Latin

# PGLS with caper package ----
# Results in 

bfly <- comparative.data(tr, trends, 
                         names.col = Latin,
                         vcv=TRUE, na.omit=FALSE, warn.dropped=TRUE)
bfly.pgls<-pgls(trend ~ 1, data=bfly, lambda="ML")
summary(bfly.pgls)
bfly.lm<-gls(trend ~ 1, data=trends)


pglslist <- list()
pglslist[[1]] <- try(pgls(trend ~ local_adapt, data=bfly, lambda="ML"))
pglslist[[2]] <- try(pgls(trend ~ annual_plast, data=bfly, lambda="ML"))
pglslist[[3]] <- try(pgls(trend ~ lg_sim, data=bfly, lambda="ML"))
pglslist[[4]] <- try(pgls(trend ~ zordspec, data=bfly, lambda="ML"))
pglslist[[5]] <- try(pgls(trend ~ zlastspec, data=bfly, lambda="ML"))
pglslist[[6]] <- try(pgls(trend ~ LG_trend, data=bfly, lambda="ML"))
pglslist[[7]] <- try(pgls(trend ~ LG_trend + lg_sim, data=bfly, lambda="ML"))
pglslist[[8]] <- try(pgls(trend ~ LG_trend * lg_sim, data=bfly, lambda="ML"))

names(pglslist) <- c("local_adapt", "annual_plast", "lg_sim", "zordspec", "zlastspec", "LG_trend", "additive", "interaction")

pglslist <- pglslist[lapply(pglslist, class) == "pgls"]

coefs_noboot <- pglslist %>% 
  purrr::map(.f = function(x){
    df <- data.frame(summary(x)$coefficients) %>% 
      mutate(var = dimnames(summary(x)$coefficients)[[1]],
             lambda = summary(x)$param[2],
             adjR2 = summary(x)$adj.r.squared)
    return(df)
  }) %>% 
  bind_rows(.id = "id")

write.csv(coefs_noboot, "data/pgls_noboot.csv", row.names = FALSE)


# phylogenetic bootstrap ----
out <- list()
for (i in 1:1000){
  tmpdat <- data.frame(cordat[,,i])
  names(tmpdat) <- names(traits)[c(8,9,2,3,4,6,5)]
  
  trends <- tmpdat %>% mutate(CommonName = sp)
  
  best.tree <- read.nexus(file = "data/RAxML_bipartitions.OHbflyBOLD_JRA_cds_aligned.gBlocks.consensus.famsConstrained.nex")
  tr <- root(best.tree, outgroup = 49, resolve.root = TRUE) #Bombyxmori outgroup, silk moth
  tr$tip.label <- gsub(pattern = "'", replacement = "", x = tr$tip.label, fixed = TRUE)
  
  species <- read.csv("data/species_names.csv") %>% 
    dplyr::select(CommonName, Genus, Species) %>%
    distinct() %>% 
    mutate(Latin = paste(Genus, Species, sep = " ")) 
  trends <- left_join(trends, species)
  
  spp <- read.csv("data/OHbflyID.csv")
  trends <- trends %>% 
    left_join(distinct(spp[,c("CommonName", "species_name")])) %>% 
    mutate(species_name = gsub(pattern = " ", replacement = "", x = species_name, fixed = TRUE))
  trends$species_name[which(trends$CommonName %in% c("Azures", "Spring/Summer Azure"))] <- "Celastrinaladon"
  
  trends <- trends[which(trends$species_name %in% tr$tip.label), ]
  trends <- droplevels(trends)
  newtips <- trends$Latin[match(tr$tip.label, trends$species_name)]
  tr$tip.label <- newtips
  
  trend <- trends %>% pull(trend)
  names(trend) <- trends %>% pull(Latin)
  tr <- keep.tip(tr, tip = names(trend))
  
  trends <- data.frame(trends)
  row.names(trends) <- trends$Latin
  
  # PGLS with caper package ----
  # Results in Table S3
  
  bfly <- comparative.data(tr, trends, 
                           names.col = Latin,
                           vcv=TRUE, na.omit=FALSE, warn.dropped=TRUE)
  
  pglslist <- list()
  pglslist[[1]] <- try(pgls(trend ~ local_adapt, data=bfly, lambda="ML"))
  pglslist[[2]] <- try(pgls(trend ~ annual_plast, data=bfly, lambda="ML"))
  pglslist[[3]] <- try(pgls(trend ~ lg_sim, data=bfly, lambda="ML"))
  pglslist[[4]] <- try(pgls(trend ~ zordspec, data=bfly, lambda="ML"))
  pglslist[[5]] <- try(pgls(trend ~ zlastspec, data=bfly, lambda="ML"))
  pglslist[[6]] <- try(pgls(trend ~ LG_trend, data=bfly, lambda="ML"))
  pglslist[[7]] <- try(pgls(trend ~ LG_trend + lg_sim, data=bfly, lambda="ML"))
  pglslist[[8]] <- try(pgls(trend ~ LG_trend * lg_sim, data=bfly, lambda="ML"))
  
  names(pglslist) <- c("local_adapt", "annual_plast", "lg_sim", "zordspec", "zlastspec", "LG_trend", "additive", "interaction")
  
  pglslist <- pglslist[lapply(pglslist, class) == "pgls"]
  
  coefs <- pglslist %>% 
    purrr::map(.f = function(x){
      df <- data.frame(summary(x)$coefficients) %>% 
        mutate(var = dimnames(summary(x)$coefficients)[[1]],
               lambda = summary(x)$param[2],
               adjR2 = summary(x)$adj.r.squared)
      return(df)
    }) %>% 
    bind_rows(.id = "id")
  
  out[[i]] <- coefs %>% mutate(sim = i)
  
}
bsout <- bind_rows(out)

saveRDS(bsout, "pgls_bootstrap.rds")

bsout <- readRDS("data/pgls_bootstrap.rds")

coefs_boot <- bsout %>% group_by(id, var) %>% summarise(med = quantile(Estimate, .5),
                                          lo = quantile(Estimate, .025),
                                          hi = quantile(Estimate, .975),
                                          t_med = quantile(t.value, .5),
                                          t_lo = quantile(t.value, .025),
                                          t_hi = quantile(t.value, .975),
                                          p_med = quantile(Pr...t.., .5),
                                          p_lo = quantile(Pr...t.., .025),
                                          p_hi = quantile(Pr...t.., .975))
write.csv(coefs_boot, "data/pgls_boot.csv", row.names = FALSE)



species <- read.csv("data/species_names.csv") %>% 
  dplyr::select(CommonName, Genus, Species) %>%
  distinct() %>% 
  mutate(Latin = paste(Genus, Species, sep = " ")) 
plttraits <- left_join(traits, species)

# Fig. S4 ----

ggplot(plttraits, aes(x = LG_trend, y = trend, label = Latin)) +
  geom_point() +
  geom_text_repel(size = 4, fontface = 'italic') +
  # scale_y_continuous(limits=c(0,100), expand = expansion(mult = c(0, 0))) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .25) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .25) +
  geom_abline(intercept = -0.021, slope = 0.502, linetype = "solid", alpha = .4) +
  annotate("text", x = 0.04, y = .1, label = "paste(italic(t), \" = 2.84, \", d.f., \" = 28, \",italic(P), \" = 0.0083\")", parse = TRUE) +
  annotate("text", x = 0.04, y = .09, label = "paste(italic(R)^{2}, \" = 0.20\")", parse = TRUE) +
  xlab("Observed trend in last generation size") +
  ylab("Observed statewide population trend\n(1st generation, 1996-2022)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "figS4.tif", path = "figures", device='tiff', dpi=600)

# Fig. S5 ----

ggplot(plttraits, aes(x = lg_sim, y = trend, label = Latin)) +
  geom_point() +
  geom_text_repel(size = 4, fontface = 'italic') +
  # scale_y_continuous(limits=c(0,100), expand = expansion(mult = c(0, 0))) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .25) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .25) +
  geom_abline(intercept = -0.019, slope = 0.0548, alpha = .4, linetype = "solid") +
  annotate("text", x = -.35, y = .1, label = "paste(italic(t), \" = 2.15, \", d.f., \" = 28, \",italic(P), \" = 0.0396\")", parse = TRUE) +
  annotate("text", x = -.35, y = .09, label = "paste(italic(R)^{2}, \" = 0.11\")", parse = TRUE) +
  xlab("Simulated change in overwinter population growth\nwith 1 std. dev. larger last generation") +
  ylab("Observed statewide population trend\n(1st generation, 1996-2022)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = "figS5.tif", path = "figures", device='tiff', dpi=600)


