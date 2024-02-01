# Cool phenology figure

source('code/01_data_prep.R')

theme_set(theme_bw(base_size = 14)) 

genpop <- readRDS("data/genpops.rds")
moddat <- readRDS("data/modeling_data.rds")

# Filter data, add last generation size variables ----
popdat <- genpop %>% 
  filter(CommonName %in% moddat$CommonName) %>% 
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


quantile_df <- function(x, probs = c(0.33, 0.66)) {
  tibble(
    val = quantile(x, probs, na.rm = TRUE),
    quant = probs
  )
}

genprop <- popdat %>% 
  group_by(CommonName) %>% 
  reframe(quantile_df(lastprop)) %>% 
  mutate(volt = rep(c("low", "high"), 30)) %>% 
  pivot_wider(id_cols = CommonName, names_from = volt, values_from = val)

voltprop <- popdat %>% 
  left_join(genprop) %>% 
  ungroup() %>% 
  mutate(volt = case_when(lastprop < low ~ "low",
                          lastprop > high ~ "high",
                          TRUE ~ "mid")) %>% 
  group_by(CommonName, gen, volt) %>% 
  summarise(prop_gen = sum(trap.N, na.rm = TRUE) / sum(trap.N.allgens, na.rm = TRUE),
            mu_gen_doy = weighted.mean(trap.mu.doy, trap.N.allgens, na.rm = TRUE),
            mu_gen_gdd = weighted.mean(trap.mu.gdd, trap.N.allgens, na.rm = TRUE),
            mu_date = as.Date(as.Date("2019-12-31",  "%Y-%m-%d") + ddays(mu_gen_doy)))


species <- read.csv("data/species_names.csv") %>% 
  dplyr::select(CommonName, Genus, Species) %>%
  distinct() %>% 
  mutate(Latin = paste(Genus, Species, sep = " ")) 

voltprop <- left_join(voltprop, species) %>% 
  group_by(CommonName) %>% 
  mutate(sp_mean = mu_gen_doy[which(gen == 1 & volt == "mid")])

sp_levels <- voltprop %>% arrange(sp_mean) %>% pull(Latin) %>% unique()

voltprop$Latin_ordered = factor(voltprop$Latin, levels=sp_levels, ordered = TRUE)

    
ggplot(voltprop %>% filter(volt != "mid"), aes(x = Latin_ordered, y = mu_date, group = volt, color = gen)) +
  geom_point(aes(size = prop_gen), shape = 15, position = position_dodge(width = 1)) +
  scale_size_area(name = "Proportion:", breaks = c(0.25, 0.75)) +
  coord_flip() +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_brewer(name = "Generation:", palette = "Dark2") +
  geom_vline(xintercept = seq(1.5, 29.5,length.out = 29), alpha = .5) +
  geom_vline(xintercept = seq(1, 30,length.out = 30), alpha = .5, linetype = "dotted") +
  ylab("Peak phenology date of generation") +
  xlab(NULL) +
  ggtitle("Voltinism variation by Ohio butterflies", subtitle = "Rows by upper/lower tercile of observed last generation size") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.text.y = element_text(face = "italic")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(plot.subtitle = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom") +
  theme(
        legend.margin = margin(0, 5, 5, 5),
        legend.spacing.x = unit(1, "mm"),
        legend.spacing.y = unit(1, "mm")) +
  guides(color = guide_legend(override.aes = list(size=5)))

ggsave(filename = "figS2.tif", path = "figures", device='tiff', dpi=600, width = 7.5, height = 13, units = "in")

  

