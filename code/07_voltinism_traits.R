## Header ---------------------------
## Script name: 07_voltinism_traits.R
## Purpose of script: 
## Author: Tyson Wepprich
## Date Created: 2024-01-22
## License: CC0 1.0 Universal
## Email: tyson.wepprich@gmail.com
## ---
## Notes: 
## 
## ---


source('01_data_prep.R')

# Methods verification
# shows how GDD always had lower class uncertainty that DOY
mmcompare <- readRDS("mixmodcomparison.rds")
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




allpars <- readRDS("species_models.rds")

# simulated benefit of last generation
# 8 species not enough data unless filter changed from 10/site to 5/site
# 8 negative, only EuroSkip & Lil Glassywing signif.
# 22 positive, 14 signif
lg_sim1 <- allpars %>% filter(model == "lg_sim1")
lg_sim1 %>% filter(!is.na(std.error)) %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)

# simulated slope of mean site LG on LG benefit to OW lambda
# 9 species significant, only 1 with negative and it doesn't have zlastann in its selected OW model
lg_sim2 <- allpars %>% filter(model == "lg_sim2") %>% 
  filter(p.value < 0.05 & term == "zlastsite")


zlastsite_slope <- allpars %>% filter(model == "lg_sim2") %>% 
  pivot_wider(id_cols = CommonName, names_from = term, values_from = estimate) %>% 
  arrange(zlastsite)

ggplot(zlastsite_slope) + geom_abline(aes(intercept = `(Intercept)`, slope = zlastsite), alpha = .5) +
  scale_x_continuous(limits = c(-5, 5)) +
  scale_y_continuous(limits = c(-5, 5))


# local adaptation
allpars %>% filter(model == "vdp_2") %>% 
  filter(term == "zordsite") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)
# pheno plasticity
allpars %>% filter(model == "vdp_2") %>% 
  filter(term == "zordann") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)
# signif co-gradient
sp1 <- allpars %>% filter(model == "vdp_2") %>% 
  filter(term == "zordann") %>% filter(p.value < 0.05) %>%  arrange(estimate) %>%  data.frame() %>% round_df(digits = 2) %>%  pull(CommonName)
sp2 <- allpars %>% filter(model == "vdp_2") %>% 
  filter(term == "zordsite") %>% filter(p.value < 0.05) %>%  arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)%>%  pull(CommonName)

# ow model
allpars %>% filter(model == "ow") %>% 
  filter(term == "zdensann") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)

# trend
allpars %>% filter(model == "trend") %>% 
  filter(term == "zyear") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 3)



# compare some traits with trends
traits <- allpars %>% 
  group_by(CommonName) %>% 
  summarise(local_adapt = estimate[which(term == "zordsite" & model == "vdp_2")],
            annual_plast = estimate[which(term == "zordann" & model == "vdp_2")],
            trend = estimate[which(term == "zyear" & model == "trend")],
            lg_sim = ifelse(is.na(std.error[which(term == "(Intercept)" & model == "lg_sim1")]), NA, estimate[which(term == "(Intercept)" & model == "lg_sim1")]))
traits <- left_join(traits, distinct(moddat2[,c("CommonName", "maxbrood", "zordspec", "zlastspec")]))

# trends <- readRDS("trendbygen.rds")
# traits <- left_join(traits, trends[which(trends$gen == "allbut"),c("CommonName", "trend")])
# GGally::ggpairs(traits[, -1])

# get trend in LG for all species without backwards variable selection results 
outlist <- list()
sp <- sort(unique(moddat2$CommonName))
for (i in sp){
  print(i)
  tmp <- moddat2 %>% filter(CommonName == i) %>% 
    droplevels()
  
  mod <- lmer(lastloglambda ~ zdensann + zyear + (zordsite + zordann)^2 + 
                (1|SiteID),
              data = tmp, 
              control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e6)))
  
  # step_fm <- step(mod, alpha.fixed = .1)
  # mod <- get_model(step_fm)
  print(summary(mod))
  # print(try(ggpredict(mod, c("zyear", "zordann")) %>% plot() + ggtitle(i)))
  # print(try(ggpredict(mod, c("zordsite", "zordann")) %>% plot() + ggtitle(i)))
  if (class(mod) == "lmerModLmerTest"){
    r2 <- MuMIn::r.squaredGLMM(mod)
    lastgen_param <- tidy(mod) %>% 
      mutate(model = "lastgen",
             # overdisp = blmeco::dispersion_glmer(mod),
             marg_r2 = r2[1,1],
             cond_r2 = r2[1,2])
  }else{
    r2 <- glance(mod)$r.squared
    lastgen_param <- tidy(mod) %>% 
      mutate(model = "lastgen",
             # overdisp = NA,
             marg_r2 = r2,
             cond_r2 = NA)
  }
  outlist[[(length(outlist)+1)]] <- lastgen_param %>% mutate(CommonName = i)
  
}

lgtrend <- bind_rows(outlist) %>% filter(term == "zyear") %>% 
  select(estimate, CommonName) %>% 
  mutate(LG_trend = estimate) %>% 
  select(-estimate)

traits <- left_join(traits, lgtrend)

bind_rows(outlist) %>% filter(term == "zyear") %>% arrange(estimate) %>%  data.frame() %>% round_df(digits = 2)
# traits <- traits %>% filter(CommonName != "Gemmed Satyr") #outlier test

medtrend = median(traits$trend)
medsim = median(traits$lg_sim)
medlgtrend = median(traits$LG_trend)

traits %>% filter(lg_sim < medsim, LG_trend > medlgtrend) # LOST GENERATIONS? median or zero?




# species trait table ----
species <- read.csv("data/species_names.csv") %>% 
  dplyr::select(CommonName, Genus, Species) %>%
  distinct() %>% 
  mutate(Latin = paste(Genus, Species, sep = " ")) 
species_table <- left_join(species_table, species)

# LG Trend
lgtrend_trait <- bind_rows(outlist) %>% filter(term == "zyear") %>% 
  select(estimate, std.error, CommonName) %>% 
  mutate(LG_trend = paste0(round(estimate,3), "  [", round(estimate-1.96*std.error,3), ", ", round(estimate+1.96*std.error,3), "]"))
species_table <- left_join(species_table, lgtrend_trait[,c("CommonName", "LG_trend")])

modtraits <- allpars %>% 
  group_by(CommonName) %>% 
  summarise(local_adapt = estimate[which(term == "zordsite" & model == "vdp_2")],
            local_adapt_se = std.error[which(term == "zordsite" & model == "vdp_2")],
            annual_plast = estimate[which(term == "zordann" & model == "vdp_2")],
            annual_plast_se = std.error[which(term == "zordann" & model == "vdp_2")],
            trend = estimate[which(term == "zyear" & model == "trend")],
            trend_se = std.error[which(term == "zyear" & model == "trend")],
            lg_sim = estimate[which(term == "(Intercept)" & model == "lg_sim1")],
            lg_sim_se = std.error[which(term == "(Intercept)" & model == "lg_sim1")]) %>% 
  mutate(LGspatial = paste0(round(local_adapt,3), "  [", round(local_adapt-1.96*local_adapt_se,3), ", ", round(local_adapt+1.96*local_adapt_se,3), "]"),
         LGannual = paste0(round(annual_plast,3), "  [", round(annual_plast-1.96*annual_plast_se,3), ", ", round(annual_plast+1.96*annual_plast_se,3), "]"),
         LGsim = paste0(round(lg_sim,3), "  [", round(lg_sim-1.96*lg_sim_se,3), ", ", round(lg_sim+1.96*lg_sim_se,3), "]"),
         PopTrend = paste0(round(trend,3), "  [", round(trend-1.96*trend_se,3), ", ", round(trend+1.96*trend_se,3), "]"))
species_table <- left_join(species_table, modtraits[,c("CommonName", "LGspatial", "LGannual", "LGsim", "PopTrend")])
species_table$Sample <- paste(species_table$uniqSY, species_table$uniqSYlam, sep = "/")

species_table_out <- species_table[which(!is.na(species_table$LG_trend)), c("Latin", "Sample", "Generations", "meandoy", "LG_trend", "LGspatial", "LGannual", "LGsim", "PopTrend")] %>% 
  arrange(Latin)

write.csv(species_table_out, "Table1.csv", row.names = FALSE)

traits_clean <- traits
names(traits_clean) <- c("CommonName", "LG response to\nsite phenology", "LG response to\nannual phenology", "Population\ntrend", 
                         "Simulated effect\nof larger LG", "Max # of\ngenerations", "Mean penult.\npeak date", 
                         "Mean\nLG size", "Trend in\nLG size")
traits_clean <- traits_clean[,c(1,6,7,8, 9, 2,3,5,4)]
GGally::ggpairs(traits_clean[, -1])

M <- cor(as.matrix(traits_clean[,-c(1:2)]))
testRes = cor.mtest(M, conf.level = 0.95)
corrplot(M, type = "upper", addCoef.col = 'black')
# 
species <- read.csv("data/species_names.csv") %>% 
  dplyr::select(CommonName, Genus, Species) %>%
  distinct() %>% 
  mutate(Latin = paste(Genus, Species, sep = " ")) 
plttraits <- left_join(traits, species)

ggplot(plttraits, aes(x = lg_sim, y = trend, label = Latin)) +
  geom_point() +
  geom_text_repel(size = 4, fontface = 'italic') +
  # scale_y_continuous(limits=c(0,100), expand = expansion(mult = c(0, 0))) +
  # geom_hline(yintercept = meantrend, linetype = "dashed", alpha = .4) +
  # geom_vline(xintercept = meansim, linetype = "dashed", alpha = .4) +
  geom_abline(intercept = -0.0181, slope = 0.0513, linetype = "dashed", alpha = .3) +
  annotate("text", x = -.4, y = .1, label = "paste(italic(t), \" = 2.03, \", d.f., \" = 28, \",italic(P), \" = 0.052\")", parse = TRUE) +
  annotate("text", x = -.4, y = .09, label = "paste(italic(R)^{2}, \" = 0.097\")", parse = TRUE) +
  xlab("Simulated change in overwinter population growth\nwith 1 std. dev. larger last generation") +
  ylab("Observed statewide population trend\n(1st generation, 1996-2022)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(plttraits, aes(x = LG_trend, y = trend, label = Latin)) +
  geom_point() +
  geom_text_repel(size = 4, fontface = 'italic') +
  # scale_y_continuous(limits=c(0,100), expand = expansion(mult = c(0, 0))) +
  # geom_hline(yintercept = meantrend, linetype = "dashed", alpha = .4) +
  # geom_vline(xintercept = meansim, linetype = "dashed", alpha = .4) +
  geom_abline(intercept = -0.0217, slope = 0.5227, linetype = "solid", alpha = .3) +
  annotate("text", x = -.01, y = .1, label = "paste(italic(t), \" = 2.93, \", d.f., \" = 28, \",italic(P), \" = 0.0067\")", parse = TRUE) +
  annotate("text", x = -.01, y = .09, label = "paste(italic(R)^{2}, \" = 0.21\")", parse = TRUE) +
  xlab("Observed trend in last generation size") +
  ylab("Observed statewide population trend\n(1st generation, 1996-2022)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Figure S5
ggplot(plttraits, aes(x = LG_trend, y = lg_sim, label = Latin, color = trend)) +
  geom_point() +
  scale_color_continuous(name = "Population\ntrend") +
  geom_text_repel(size = 4) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = .4) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = .4) +
  xlab("Observed trend in last generation size") +
  ylab("Simulated change in overwinter population growth\nwith 1 std. dev. larger last generation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(traits, aes(x = LG_trend, y = trend, label = CommonName, color = lg_sim)) +
  geom_point() +
  geom_text_repel(size = 4) +
  xlab("Observed trend in last generation size") +
  ylab("Simulated change in overwinter population growth\nwith 1 std. dev. larger last generation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# correlations with bootstrap----
set.seed(123)

# from above model results with LG trend forced to fit
lgtrend <- bind_rows(outlist) %>% filter(term == "zyear") 

cordat <- array(data = NA, dim = c(30, 7, 1000))
for (i in 1:30){
  
  tmpars <- allpars %>% filter(CommonName == sp[i])
  
  # species phenology
  cordat[i,1,] <- sur::boot.mean(moddat2$zord[moddat2$CommonName == sp[i]], 1000)$bootstrap.samples
  # species LG size
  cordat[i,2,] <- sur::boot.mean(moddat2$lastloglambda[moddat2$CommonName == sp[i]], 1000)$bootstrap.samples
  # species LG trend
  cordat[i,3,] <- rnorm(n = 1000, mean = lgtrend$estimate[lgtrend$CommonName == sp[i]], sd = lgtrend$std.error[lgtrend$CommonName == sp[i]])
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
corrplot::corrplot(Mcoef, type = "upper", addCoef.col = 'black', tl.col	= "black", col = COL2('BrBG', 100))

conf <- paste0("(", format(Mlo, digits=1), ",", format(Mhi, digits=1), ")")
xs <- row(Mlo)
ys <- (ncol(Mlo)+1) - col(Mlo)
text(xs[lower.tri(xs)], ys[lower.tri(xs)], conf[lower.tri(xs)], pos=1, cex=.75)


trait.bs <- apply(cordat, c(1,2), FUN = stats::quantile, probs = c(0.025, .5, 0.975))
trait.bs.lo <- data.frame(CommonName = sp, trait.bs[1,,], quant = "lo")
trait.bs.med <- data.frame(CommonName = sp, trait.bs[2,,], quant = "med")
trait.bs.hi <- data.frame(CommonName = sp, trait.bs[3,,], quant = "hi")
trait.df <- rbind(trait.bs.lo, trait.bs.med, trait.bs.hi)
names(trait.df) <- c(names(traits)[c(1,7,8,9,2,3,5,4)], "quant")

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

# bootstrap intriguing interaction
out <- list()
# cordat <- cordat[-c(6, 24, 9, 4, 5, 10),,]
for (i in 1:1000){
  tmpdat <- data.frame(cordat[,,i])
  names(tmpdat) <- names(traits)[c(7,8,9,2,3,5,4)]
  mod <- lm(trend ~ LG_trend * lg_sim, data = tmpdat)
  # pmod <- pgls(trend ~ LG_trend*lg_sim - 1, data=bfly, lambda="ML")
  out[[i]] <- tidy(mod) %>% mutate(sim = i)
}
bsout <- bind_rows(out)
bsout %>% group_by(term) %>% summarise(med = quantile(estimate, .5),
                                       lo = quantile(estimate, .025),
                                       hi = quantile(estimate, .975),
                                       p_med = quantile(p.value, .5),
                                       p_lo = quantile(p.value, .025),
                                       p_hi = quantile(p.value, .975))

# bootstrap intriguing interaction with phylo
out <- list()


# cordat <- cordat[-c(6, 24, 9, 4, 5, 10),,]
for (i in 1:1000){
  tmpdat <- data.frame(cordat[,,i])
  names(tmpdat) <- names(traits)[c(7,8,9,2,3,5,4)]
  # mod <- lm(trend ~ LG_trend * lg_sim, data = tmpdat)
  
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
  # Results in S1 Appendix Tables A and B
  
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
  
  # out[[i]] <- tidy(mod) %>% mutate(sim = i)
}
bsout <- bind_rows(out)

# bsout %>% group_by(term) %>% summarise(med = quantile(estimate, .5),
#                                        lo = quantile(estimate, .025),
#                                        hi = quantile(estimate, .975),
#                                        p_med = quantile(p.value, .5),
#                                        p_lo = quantile(p.value, .025),
#                                        p_hi = quantile(p.value, .975))


saveRDS(bsout, "pgls_bootstrap.rds")

bsout <- readRDS("pgls_bootstrap.rds")

bsout %>% group_by(id, var) %>% summarise(med = quantile(Estimate, .5),
                                          lo = quantile(Estimate, .025),
                                          hi = quantile(Estimate, .975),
                                          t_med = quantile(t.value, .5),
                                          t_lo = quantile(t.value, .025),
                                          t_hi = quantile(t.value, .975),
                                          p_med = quantile(Pr...t.., .5),
                                          p_lo = quantile(Pr...t.., .025),
                                          p_hi = quantile(Pr...t.., .975))


traits_cut <- traits[-c(6, 24, 9, 4, 5, 10),]
lmfit <- lm(trend ~ lg_sim, data = traits)
lmfit <- lm(trend ~ LG_trend, data = traits)
lmfit <- lm(trend ~ zlastspec, data = traits)
lmfit <- lm(trend ~ LG_trend + lg_sim, data = traits)
lmfit <- lm(trend ~ LG_trend*lg_sim, data = traits)
# lmfit1 <- lm(trend ~ LG_trend + lg_sim, data = traits_cut)
# lmfit2 <- lm(trend ~ LG_trend*lg_sim, data = traits_cut)


ggpredict(lmfit, c("lg_sim[-0.6, .8]", "LG_trend[-0.01,0.05]")) %>% plot(add.data = TRUE)

pltdat <- ggpredict(lmfit, c("lg_sim[-0.6, .7]", "LG_trend[-0.0135,0.0565]")) %>% plot(add.data = TRUE)
pltdat <- ggpredict(lmfit, c("lg_sim[-0.6, .6]", "LG_trend[-0.0135,0.0565]")) %>% plot(add.data = TRUE, color = "eight")

# change attributes so facets labeled (3rd one only?)
attr(pltdat$data, "terms") <- c("Simulated effect of larger LG", "Last generation trend")

pltdat + 
  # scale_colour_brewer(palette = "Set1", direction = -1, labels = c("-1 SD", "+1 SD")) +
  # scale_fill_brewer(palette = "Set1", direction = -1) +
  # scale_y_continuous(limits=c(-2,1), expand = expansion(mult = c(0, 0))) +
  labs(title = NULL, x = "Overwinter effect of larger last generation", y = "Statewide population trend") +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
# guides(color=guide_legend(title="Mean site\nlast generation\nsize"))


data.frame(CommonName=sp, stdres = as.numeric(residuals(lmfit)/sqrt(var(residuals(lmfit)))))















# Population trends vs phylogeny


library(ape)
# library(phangorn)
# library(phytools)
library(caper)

# Read phylogenetic tree ----
# used in PGLS and supplementary figure

best.tree <- read.nexus(file = "data/RAxML_bipartitions.OHbflyBOLD_JRA_cds_aligned.gBlocks.consensus.famsConstrained.nex")
tr <- root(best.tree, outgroup = 49, resolve.root = TRUE) #Bombyxmori outgroup, silk moth
tr$tip.label <- gsub(pattern = "'", replacement = "", x = tr$tip.label, fixed = TRUE)
trends <- traits

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


# choose whether to include migrant species or not
# trends <- trends %>% filter(ResStatus != "immigrant")

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
# Results in S1 Appendix Tables A and B

bfly <- comparative.data(tr, trends, 
                         names.col = Latin,
                         vcv=TRUE, na.omit=FALSE, warn.dropped=TRUE)


bfly.pgls<-pgls(trend ~ 1, data=bfly, lambda="ML")
summary(bfly.pgls)

# Univariate models
pglslist <- list()
pglslist[[1]] <- pgls(trend ~ local_adapt, data=bfly, lambda="ML")
pglslist[[2]] <- pgls(trend ~ annual_plast, data=bfly, lambda="ML")
pglslist[[3]] <- pgls(trend ~ lg_sim, data=bfly, lambda="ML")
pglslist[[4]] <- pgls(trend ~ maxbrood, data=bfly, lambda="ML")
pglslist[[5]] <- pgls(trend ~ zordspec, data=bfly, lambda="ML")
pglslist[[6]] <- pgls(trend ~ zlastspec, data=bfly, lambda="ML")
pglslist[[7]] <- pgls(trend ~ LG_trend, data=bfly, lambda="ML")
pglslist[[8]] <- pgls(trend ~ LG_trend + lg_sim, data=bfly, lambda="ML")
pglslist[[9]] <- pgls(trend ~ LG_trend * lg_sim, data=bfly, lambda="ML")

AIC(bfly.pgls, pglslist[[7]], pglslist[[8]], pglslist[[9]])

coefs <- pglslist %>% 
  purrr::map(.f = function(x){
    df <- data.frame(summary(x)$coefficients) %>% 
      mutate(var = dimnames(summary(x)$coefficients)[[1]],
             lambda = summary(x)$param[2],
             adjR2 = summary(x)$adj.r.squared)
    return(df)
  }) %>% 
  bind_rows(.id = "id")

write.csv(round_df(coefs, 3), "pglstraits.csv", row.names = FALSE)

res<- residuals(pglslist[[10]], phylo = TRUE)
res<- res/sqrt(var(res))[1]
rownames(trends)[(abs(res)>2)]


profile_lambda=pgls.profile(pglslist[[9]], which="lambda") # vary lambda
plot(profile_lambda)


plot(trend ~ LG_trend, data = trends)
abline(pglslist[[7]])



# Figure 5 ----
# figure with multiple panels
# sorry, this has manual changes, each trait model selected, df changed, plot saved in list

levels(trends$ResStatus) <- c("Migratory", "Year-round")
trends$Range_Num <- factor(as.character(trends$Range_Num), levels = levels(trends$Range_Num)[c(3, 1, 2)])
levels(trends$Range_Num) <- c("Southern", "Core", "Northern")
levels(trends$HostCategory) <- c("Forb", "Graminoid", "Woody")
levels(trends$Special) <- c("Generalist", "Family", "Genus")
levels(trends$WinterStage) <- c("Adult", "Egg", "Larva", "Pupa")




mod <- pglslist[[4]]

df <- data.frame(ResStatus = levels(trends$WinterStage),
                 estimate = summary(mod)$coefficients[,1],
                 se = summary(mod)$coefficients[,2])

theme_set(theme_bw(base_size = 16)) 
# pltlist <- list()

# Discrete traits
plt <- ggplot(data = trends, aes(x = WinterStage, y = estimate)) +
  # geom_boxplot() +
  geom_jitter(alpha = .5, width = .2, height = 0) +
  geom_point(data = df, aes(x = ResStatus, y = estimate), size = 4, alpha = 1, shape = 15) +
  geom_linerange(data = df, aes(x = ResStatus, ymin = estimate - 1.96 * se, ymax = estimate + 1.96 * se), size = 1) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  xlab("Overwinter Stage") +
  ylab("Statewide abundance trend") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5))


plt
pltlist[[4]] <- plt

# ggsave(filename = "HostCategory_pgls.png", device = "png", plot = plt, width = 6, height = 6, units = "in")

p <- ggarrange(pltlist[[1]], pltlist[[2]] + rremove("ylab"), pltlist[[3]] + rremove("ylab"), 
               pltlist[[4]], pltlist[[5]]+ rremove("ylab"), pltlist[[6]]+ rremove("ylab"),
               labels = c("A", "B", "C", "D", "E", "F"),
               ncol = 3, nrow = 2, align = "v")
p
ggsave(plot = p, filename =  "traits.png", device = "png", width = 12, height = 8, units = "in")

ggsave(filename = "traits.tiff", plot = p, device = "tiff", dpi = "retina", 
       width = 12, height = 8, units = "in", compression = "lzw")


