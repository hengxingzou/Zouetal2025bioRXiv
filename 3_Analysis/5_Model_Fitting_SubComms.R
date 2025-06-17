# MODEL FITTING FOR TOP CONTRIBUTORS ONLY AND THE REST OF THE SPECIES, AND THEIR COMPARISONS
# Hengxing Zou

########## Initialization ##########


# Double check to make sure the directory is correct

source("3_Analysis/2_Model_Fitting_Base.R")

# Change the directory according to your working environment

generated_dir = "GeneratedData/"
database_dir = "Database/"
figure_dir = "Figures/"
stats_dir = "Stats/"


########## Read Subcommunities Data ##########


# The following data are generated from 4_Analyze_Contribution.R
# Read community-weighted metrics of species with the largest contributions

Weighted_TopSpp = read_csv(paste0(generated_dir, "CWMetrics_TopSpp.csv")) %>% 
  rename("Hand.Wing.Index" = "Hand-Wing.Index")

# Read community-weighted metrics of species with the largest contributions

Weighted_RestSpp = read_csv(paste0(generated_dir, "CWMetrics_RestSpp.csv")) %>% 
  rename("Hand.Wing.Index" = "Hand-Wing.Index")

# Join CW Metrics with climate and anthrome data

cwm_topspp = Weighted_TopSpp %>% 
  inner_join(Bioclim, by = c("Region" = "ST_12", "year" = "year")) %>% 
  inner_join(Anthromes, by = c("Region" = "ST_12", "year" = "year"))

cwm_mean_topspp = cwm_topspp %>% 
  filter(Metric == "CWM")

cwm_restspp = Weighted_RestSpp %>% 
  inner_join(Bioclim, by = c("Region" = "ST_12", "year" = "year")) %>% 
  inner_join(Anthromes, by = c("Region" = "ST_12", "year" = "year"))

cwm_mean_restspp = cwm_restspp %>% 
  filter(Metric == "CWM")


########## Transform Data ##########


# Scale and transform data, top contributors

cwm_mean_topspp_t = cwm_mean_topspp %>% 
  
  # standardize years
  mutate(year_since_1969 = year - 1969) %>% 
  
  # calculate temperature and precipitation seasonality
  mutate(temp_seasonality = bio5 - bio6, precip_seasonality = bio13 - bio14) %>% 
  
  # remove unnecessary bioclim variables
  select(-(bio1:bio9), -(bio11:bio17), -bio19) %>% 
  
  # scale bioclim and land cover variables
  mutate(across(c(bio10, bio18, temp_seasonality, precip_seasonality, 
                  Settlements:Wild), ~c(scale(.)))) %>%
  
  # log transform some CWMs
  mutate(log_Mass = log(Mass), 
         log_clutch_size = log(litter_or_clutch_size_n), 
         log_GenLength = log(GenLength)) %>% 
  
  # scale all CWMs
  mutate(across(Mass:GenLength, ~c(scale(.)))) %>% 
  
  # scale log(body mass) for regression as a predictor variable
  mutate(Scaled_Mass = c(scale(log_Mass, scale = T, center = T))) %>% 
  
  # scale corrected generation lengths
  mutate(corr_GenLength = c(scale(corr_GenLength, scale = T, center = T)))

# Scale and transform data, rest of the communities

cwm_mean_restspp_t = cwm_mean_restspp %>% 
  
  # standardize years
  mutate(year_since_1969 = year - 1969) %>% 
  
  # calculate temperature and precipitation seasonality
  mutate(temp_seasonality = bio5 - bio6, precip_seasonality = bio13 - bio14) %>% 
  
  # remove unnecessary bioclim variables
  select(-(bio1:bio9), -(bio11:bio17), -bio19) %>% 
  
  # scale bioclim and land cover variables
  mutate(across(c(bio10, bio18, temp_seasonality, precip_seasonality, 
                  Settlements:Wild), ~c(scale(.)))) %>%
  
  # log transform some CWMs
  mutate(log_Mass = log(Mass), 
         log_clutch_size = log(litter_or_clutch_size_n), 
         log_GenLength = log(GenLength)) %>% 
  
  # scale all CWMs
  mutate(across(Mass:GenLength, ~c(scale(.)))) %>% 
  
  # scale log(body mass) for regression as a predictor variable
  mutate(Scaled_Mass = c(scale(log_Mass, scale = T, center = T))) %>% 
  
  # scale corrected generation lengths
  mutate(corr_GenLength = c(scale(corr_GenLength, scale = T, center = T)))


########## CWM Year Model ##########


# Year effect, top contributors

year_models_topspp = foreach(tr = md_traits) %dopar% {
  
  formula_full = as.formula(paste(tr, " ~ year_since_1969 + (1 | Region)"))
  fit = lmerTest::lmer(formula_full, data = cwm_mean_topspp_t, weights = Community_Abundance)
  fit
  
}

names(year_models_topspp) = md_traits

year_params_topspp = foreach(i = 1:length(year_models_topspp), .combine = rbind, .packages = "tidyverse") %dopar% {
  
  m = as.data.frame(summary(year_models_topspp[[i]])$coeff) %>%
    rownames_to_column("parameters") %>%
    mutate(trait = names(year_models_topspp)[i])
  colnames(m) = c("parameters", "estimate", "se", "df", "t_value", "p_value", "trait")
  
  m
  
}

# Year effect, rest of the communities

year_models_restspp = foreach(tr = md_traits) %dopar% {
  
  formula_full = as.formula(paste(tr, " ~ year_since_1969 + (1 | Region)"))
  fit = lmerTest::lmer(formula_full, data = cwm_mean_restspp_t, weights = Community_Abundance)
  fit
  
}

names(year_models_restspp) = md_traits

year_params_restspp = foreach(i = 1:length(year_models_restspp), .combine = rbind, .packages = "tidyverse") %dopar% {
  
  m = as.data.frame(summary(year_models_restspp[[i]])$coeff) %>%
    rownames_to_column("parameters") %>%
    mutate(trait = names(year_models_restspp)[i])
  colnames(m) = c("parameters", "estimate", "se", "df", "t_value", "p_value", "trait")
  
  m
  
}

# Combine and compare all subcommunities

year_all_params = rbind(year_params_cwm %>% mutate(Community = "All"), 
                        year_params_topspp %>% mutate(Community = "TopSpp"), 
                        year_params_restspp %>% mutate(Community = "RestSpp"))

# Save stats

write_csv(year_all_params, paste0(stats_dir, "Year_All_Params.csv"))

# Visualization

p_yr_all = year_all_params %>% 
  
  # Add significance values
  mutate(signif = if_else(p_value >= 0.05, F, T), 
         positive = if_else(estimate > 0, T, F), 
         trait = recode(trait, !!!trait_labels)) %>% 
  mutate(trait = factor(trait, rev(trait_labels))) %>% 
  
  # Filter only slopes
  filter(parameters != "(Intercept)") %>% 
  
  ggplot(aes(y = trait, x = estimate, color = positive, shape = Community)) + 
  geom_pointrange(aes(xmin = estimate-se, xmax = estimate+se), 
                  position = ggstance::position_dodgev(height = 0.75),
                  size = 0.5) + 
  geom_vline(xintercept = 0, color = "black") + 
  xlab("Estimate") + 
  ylab("Trait") +
  scale_color_manual(values = color_2, guide = "none") +
  scale_shape_manual(values = c(16, 17, 15), labels = c("Base", "Rest", "Top")) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.position = "bottom"
  )

p_yr_all

# Save figure

ggsave(paste0(figure_dir, "CWM_lm_Year_Compare.png"), p_yr_all, device = "png", width = 1600, height = 1200, unit = "px")


########## Model Fitting, Top Contributors ##########


# Fit models for top contributors

cwm_models_topspp = foreach(tr = md_traits_log) %dopar% {
  
  formula_full = as.formula(paste(tr, " ~ year_since_1969 +
                            bio10 + bio18 + temp_seasonality + precip_seasonality +
                            Settlements + Agriculture + Cultured + 
                            temp_seasonality:precip_seasonality +
                            year_since_1969:bio10 + year_since_1969:bio18 +
                            year_since_1969:temp_seasonality + year_since_1969:precip_seasonality +
                            year_since_1969:Settlements + year_since_1969:Cultured +
                            year_since_1969:Agriculture + 
                            (1 | Region)"))
  
  fit = lmerTest::lmer(formula_full, data = cwm_mean_topspp_t, weights = Community_Abundance, 
                       na.action = "na.fail")
  
  fit
  
}

names(cwm_models_topspp) = md_traits_log

cwm_params_topspp = foreach(i = 1:length(cwm_models_topspp), .combine = rbind, .packages = "tidyverse") %dopar% {
  
  m = as.data.frame(summary(cwm_models_topspp[[i]])$coeff) %>% 
    rownames_to_column("parameters") %>% 
    mutate("trait" = names(cwm_models_topspp)[i])
  colnames(m) = c("parameters", "estimate", "se", "df", "t_value", "p_value", "trait")
  
  m
  
}

# Fit models for rest of the species

cwm_models_restspp = foreach(tr = md_traits_log) %dopar% {
  
  formula_full = as.formula(paste(tr, " ~ year_since_1969 +
                            bio10 + bio18 + temp_seasonality + precip_seasonality +
                            Settlements + Agriculture + Cultured + 
                            temp_seasonality:precip_seasonality +
                            year_since_1969:bio10 + year_since_1969:bio18 +
                            year_since_1969:temp_seasonality + year_since_1969:precip_seasonality +
                            year_since_1969:Settlements + year_since_1969:Cultured +
                            year_since_1969:Agriculture + 
                            (1 | Region)"))
  
  fit = lmerTest::lmer(formula_full, data = cwm_mean_restspp_t, weights = Community_Abundance, 
                       na.action = "na.fail")
  
  fit
  
}

names(cwm_models_restspp) = md_traits_log

cwm_params_restspp = foreach(i = 1:length(cwm_models_restspp), .combine = rbind, .packages = "tidyverse") %dopar% {
  
  m = as.data.frame(summary(cwm_models_restspp[[i]])$coeff) %>% 
    rownames_to_column("parameters") %>% 
    mutate("trait" = names(cwm_models_restspp)[i])
  colnames(m) = c("parameters", "estimate", "se", "df", "t_value", "p_value", "trait")
  
  m
  
}


########## Compare Full and Subcommunities ##########


cwm_all_params = rbind(cwm_params %>% mutate(Community = "All"), 
                   cwm_params_topspp %>% mutate(Community = "TopSpp"), 
                   cwm_params_restspp %>% mutate(Community = "RestSpp"))

cwm_params_nointeractions = cwm_all_params %>% 
  
  # Add significance values
  mutate(signif = if_else(p_value >= 0.05, F, T), 
         positive = if_else(estimate > 0, T, F)) %>% 
  
  # Filter slopes
  filter(parameters != "(Intercept)") %>% 
  filter(!grepl("PC2", trait)) %>% 
  filter(!grepl(":", parameters)) %>%
  
  # Reorder variables
  mutate(parameters = recode(parameters, !!!labels)) %>% 
  mutate(parameters = factor(parameters, rev(labels))) %>%
  mutate(trait = recode(trait, !!!trait_labels)) %>% 
  mutate(trait = factor(trait, trait_labels)) %>% 
  mutate(Community = factor(Community, c("RestSpp", "TopSpp", "All")))

# Save stats

write_csv(cwm_all_params, paste0(stats_dir, "CWM_All_Params.csv"))
write_csv(cwm_params_nointeractions, paste0(stats_dir, "CWM_Params_No_Interactions.csv"))

# Visualization

p_tr_cwm_all_1 = cwm_params_nointeractions %>% 
  filter(trait %in% trait_labels[1:4]) %>% 
  ggplot(aes(y = parameters, x = estimate, color = positive, alpha = signif, shape = Community)) + 
  geom_pointrange(aes(xmin = estimate-se, xmax = estimate+se), 
                  position = ggstance::position_dodgev(height = 0.75),
                  size = 0.5) + 
  geom_vline(xintercept = 0, color = "black") + 
  scale_shape_manual(values = c(17, 15, 16)) + 
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_color_manual(values = color_2) +
  xlab("Estimate") + 
  ylab("Parameter") +
  ggh4x::facet_wrap2(. ~ trait, scales = "free_x", ncol = 4) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.position = "none"
  )

p_tr_cwm_all_2 = cwm_params_nointeractions %>% 
  filter(trait %in% trait_labels[c(8, 9, 11)]) %>% 
  ggplot(aes(y = parameters, x = estimate, color = positive, alpha = signif, shape = Community)) + 
  geom_pointrange(aes(xmin = estimate-se, xmax = estimate+se), 
                  position = ggstance::position_dodgev(height = 0.75),
                  size = 0.5) + 
  geom_vline(xintercept = 0, color = "black") + 
  scale_shape_manual(values = c(17, 15, 16)) + 
  scale_alpha_manual(values = c(0.2, 1)) +
  scale_color_manual(values = color_2) +
  xlab("Estimate") + 
  ylab("Parameter") +
  ggh4x::facet_wrap2(. ~ trait, scales = "free_x", ncol = 4) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.position = "none"
  )

p_tr_cwm_all = p_tr_cwm_all_1 / p_tr_cwm_all_2 + plot_layout(axes = "collect", design = "AAAA
                                                                                         BBB#")

p_tr_cwm_all

# Save figure

ggsave(paste0(figure_dir, "CWM_Year_Env_Compare_NoInteraction.pdf"), p_tr_cwm_all, width = 3200, height = 2400, unit = "px")

