# MODEL FITTING FOR BASE COMMUNITIES (ALL SPECIES)
# Hengxing Zou


########## Initialization ##########


# Double check to make sure the directory is correct

source("3_Analysis/1_Read_All_Data.R")

library(lmerTest)

# Figure labels for predictor variables

labels = c("year_since_1969" = "Year",
           "year_log" = "log(Year)",
           "bio10" = "Temp Warmest Q",
           "bio18" = "Precip Warmest Q",
           "temp_seasonality" = "Temp Seasonality",
           "precip_seasonality" = "Precip Seasonality",
           "Agriculture" = "Agriculture",
           "Cultured" = "Cultured Land",
           "Settlements" = "Settlements",
           "Wild" = "Wild",
           "Scaled_Mass" = "Body Mass",
           "log_Mass" = "log(Mass)",
           "temp_seasonality:precip_seasonality" = "Temp:Precip", 
           "year_since_1969:bio10" = "Year:Temp Warmest Q",
           "year_since_1969:bio18" = "Year:Precip Warmest Q",
           "year_since_1969:temp_seasonality" = "Year:Temp Sesonality",
           "year_since_1969:precip_seasonality" = "Year:Precip Sesonality",
           "year_since_1969:Agriculture" = "Year:Agriculture",
           "year_since_1969:Cultured" = "Year:Cultured Land",
           "year_since_1969:Settlements" = "Year:Settlements",
           "year_since_1969:Wild" = "Year:Wild"
)

# Change the directory according to your working environment

figure_dir = "Figures/"
stats_dir = "Stats/"


########## Transform Data ##########


# Scale and transform data, CW Mean

cwm_mean_t = cwm_mean %>% 
  
  # standardize years
  mutate(year_since_1969 = year - 1969) %>% 
  
  # calculate temperature and precipitation seasonality
  mutate(temp_seasonality = bio5 - bio6, precip_seasonality = bio13 - bio14) %>% 
  
  # remove unnecessary bioclim variables
  select(-(bio1:bio9), -(bio11:bio17), -bio19) %>% 
  
  # scale bioclim variables
  mutate(across(c(bio10, bio18, temp_seasonality, precip_seasonality), ~c(scale(.)))) %>%
  
  # scale land cover variables
  mutate(across(Settlements:Wild, ~c(scale(.)))) %>%

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


# Year effect

year_models_cwm = foreach(tr = md_traits) %dopar% {
  
  formula_full = as.formula(paste(tr, " ~ year_since_1969 + (1 | Region)"))
  
  fit = lmerTest::lmer(formula_full, data = cwm_mean_t, weights = Community_Abundance)
  
  fit
  
}

names(year_models_cwm) = md_traits

year_params_cwm = foreach(i = 1:length(year_models_cwm), .combine = rbind, .packages = "tidyverse") %dopar% {
  
  m = as.data.frame(summary(year_models_cwm[[i]])$coeff) %>%
    rownames_to_column("parameters") %>%
    mutate(trait = names(year_models_cwm)[i])
  colnames(m) = c("parameters", "estimate", "se", "df", "t_value", "p_value", "trait")
  
  m
  
}

p_yr_cwm = year_params_cwm %>% 
  
  # Add significance values
  mutate(signif = if_else(p_value >= 0.05, F, T), 
         positive = if_else(estimate > 0, T, F), 
         trait = recode(trait, !!!trait_labels)) %>% 
  mutate(trait = factor(trait, rev(trait_labels))) %>% 
  
  # Filter only slopes
  filter(parameters != "(Intercept)") %>% 
  
  ggplot(aes(y = trait, x = estimate, color = positive)) + 
  geom_pointrange(aes(xmin = estimate-se, xmax = estimate+se), 
                  size = 0.75) + 
  geom_vline(xintercept = 0, color = "black") + 
  xlab("Estimate") + 
  ylab("Trait") +
  scale_color_manual(values = color_2) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.position = "bottom"
  )

p_yr_cwm

# Save figure

ggsave(paste0(figure_dir, "CWM_lm_Year.png"), p_yr_cwm, device = "png", width = 1600, height = 1200, unit = "px")


########## CWM Year-Environment Model ##########


# Fit models

year_env_models = foreach(tr = md_traits_log) %dopar% {
  
  formula_full = as.formula(paste(tr, " ~ year_since_1969 +
                            bio10 + bio18 + temp_seasonality + precip_seasonality +
                            Settlements + Agriculture + Cultured + 
                            (1 | Region)"))
  
  fit = lmerTest::lmer(formula_full, data = cwm_mean_t, weights = Community_Abundance, 
                       na.action = "na.fail")
  
  fit
  
}

names(year_env_models) = md_traits

year_env_params = foreach(i = 1:length(year_env_models), .combine = rbind, .packages = "tidyverse") %dopar% {
  
  m = as.data.frame(summary(year_env_models[[i]])$coeff) %>% 
    rownames_to_column("parameters") %>% 
    mutate("trait" = names(year_env_models)[i])
  colnames(m) = c("parameters", "estimate", "se", "df", "t_value", "p_value", "trait")
  
  m
  
}

# Save statistics

write_csv(year_env_params, paste0(stats_dir, "/CWM_Year_Env.csv"))

# Visualization

tr_year_env = year_env_params %>% 
  
  # Add significance values
  mutate(signif = if_else(p_value >= 0.05, F, T), 
         positive = if_else(estimate > 0, T, F)) %>% 
  
  # Filter slopes
  filter(parameters != "(Intercept)") %>% 
  filter(!grepl("PC2", trait)) %>% 
  
  # Reorder variables
  mutate(parameters = recode(parameters, !!!labels)) %>% 
  mutate(parameters = factor(parameters, rev(labels))) %>%
  mutate(trait = recode(trait, !!!trait_labels)) %>% 
  mutate(trait = factor(trait, trait_labels))

p_tr_year_env_1 = tr_year_env %>% 
  filter(trait %in% trait_labels[1:4]) %>% 
  ggplot(aes(y = parameters, x = estimate, color = positive, alpha = signif)) + 
  geom_pointrange(aes(xmin = estimate-se, xmax = estimate+se), 
                  size = 0.75) + 
  geom_vline(xintercept = 0, color = "black") + 
  scale_alpha_manual(values = c(0.2, 1)) +
  xlab("Estimate") + 
  ylab("Parameter") +
  scale_color_manual(values = color_2) +
  ggh4x::facet_wrap2(. ~ trait, scales = "free_x", ncol = 4) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.position = "none"
  )

p_tr_year_env_2 = tr_year_env %>% 
  filter(trait %in% trait_labels[c(5, 6, 11)]) %>% 
  ggplot(aes(y = parameters, x = estimate, color = positive, alpha = signif)) + 
  geom_pointrange(aes(xmin = estimate-se, xmax = estimate+se), 
                  size = 0.75) + 
  geom_vline(xintercept = 0, color = "black") + 
  scale_alpha_manual(values = c(0.2, 1)) +
  xlab("Estimate") + 
  ylab("Parameter") +
  scale_color_manual(values = color_2) +
  ggh4x::facet_wrap2(. ~ trait, scales = "free_x", ncol = 4) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.position = "none"
  )

p_tr_year_env = p_tr_year_env_1 / p_tr_year_env_2 + plot_layout(axes = "collect", design = "AAAA
                                                                             BBB#")
p_tr_year_env


########## CWM Year-Environment Model ##########


# Fit models

cwm_models = foreach(tr = md_traits_log) %dopar% {
  
  formula_full = as.formula(paste(tr, " ~ year_since_1969 +
                            bio10 + bio18 + temp_seasonality + precip_seasonality +
                            Settlements + Agriculture + Cultured + 
                            temp_seasonality:precip_seasonality +
                            year_since_1969:bio10 + year_since_1969:bio18 +
                            year_since_1969:temp_seasonality + year_since_1969:precip_seasonality +
                            year_since_1969:Settlements + year_since_1969:Cultured +
                            year_since_1969:Agriculture + 
                            (1 | Region)"))
  
  fit = lmerTest::lmer(formula_full, data = cwm_mean_t, weights = Community_Abundance, 
                       na.action = "na.fail")
  
  fit
  
}

names(cwm_models) = md_traits_log

cwm_params = foreach(i = 1:length(cwm_models), .combine = rbind, .packages = "tidyverse") %dopar% {
  
  m = as.data.frame(summary(cwm_models[[i]])$coeff) %>% 
    rownames_to_column("parameters") %>% 
    mutate("trait" = names(cwm_models)[i])
  colnames(m) = c("parameters", "estimate", "se", "df", "t_value", "p_value", "trait")
  
  m
  
}

# Save statistics

write_csv(year_env_params, paste0(stats_dir, "/CWM_Full.csv"))

# Visualization

tr_cwm = cwm_params %>% 
  
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
  mutate(trait = factor(trait, trait_labels))

p_tr_cwm_1 = tr_cwm %>% 
  filter(trait %in% trait_labels[1:4]) %>% 
  ggplot(aes(y = parameters, x = estimate, color = positive, alpha = signif)) + 
  geom_pointrange(aes(xmin = estimate-se, xmax = estimate+se), 
                  size = 0.75) + 
  geom_vline(xintercept = 0, color = "black") + 
  scale_alpha_manual(values = c(0.2, 1)) +
  xlab("Estimate") + 
  ylab("Parameter") +
  scale_color_manual(values = color_2) +
  ggh4x::facet_wrap2(. ~ trait, scales = "free_x", ncol = 4) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.position = "none"
  )

p_tr_cwm_2 = tr_cwm %>% 
  filter(trait %in% trait_labels[c(8, 9, 11)]) %>% 
  ggplot(aes(y = parameters, x = estimate, color = positive, alpha = signif)) + 
  geom_pointrange(aes(xmin = estimate-se, xmax = estimate+se), 
                  size = 0.75) + 
  geom_vline(xintercept = 0, color = "black") + 
  scale_alpha_manual(values = c(0.2, 1)) +
  xlab("Estimate") + 
  ylab("Parameter") +
  scale_color_manual(values = color_2) +
  ggh4x::facet_wrap2(. ~ trait, scales = "free_x", ncol = 4) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, angle = 60, vjust = 0.5), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        strip.text = element_text(size = 10), 
        legend.position = "none"
  )

p_tr_cwm = p_tr_cwm_1 / p_tr_cwm_2 + plot_layout(axes = "collect", design = "AAAA
                                                                             BBB#")

p_tr_cwm

# Save figure

ggsave(paste0(figure_dir, "CWM_Year_Env_NoInteraction.png"), p_tr_cwm, device = "png", width = 3200, height = 2400, unit = "px")

