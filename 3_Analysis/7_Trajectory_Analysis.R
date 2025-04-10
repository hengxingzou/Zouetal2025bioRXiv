# ANALYZE TEMPORAL SHIFTS OF COMMUNITIES IN THE FUNCTIONAL SPACE
# Hengxing Zou


########## Initialization ##########


library(tidyverse)
library(patchwork)

library(ggfortify)
library(sf)
library(spdep)

library(foreach)
library(doSNOW)

# None of the trajectory analyses that use this package made it into the manuscript, 
# but we are including it just for the completeness of the code

library(trajr)

cl = makeCluster(4)
registerDoSNOW(cl)

# Double check to make sure the directory is correct

source("3_Analysis/calculate_traj.R")
source("3_Analysis/calculate_traj_metrics.R")
source("3_Analysis/plot_maps.R")
source("3_Analysis/calculate_spatial_heterogeneity.R")

# Change the directory according to your working environment

generated_dir = "GeneratedData/"
database_dir = "Database/"
figure_dir = "Figures/"
stats_dir = "Stats/"

# Read latlong categorized by BCR

category_by_bcr = read_csv(paste0(database_dir, "Grids_by_BCR.csv"))

# Read latlong map, then combine with categories

latlong_shp = read_sf(paste0(database_dir, "BBS_LatLong_strata.shp")) %>% 
  full_join(category_by_bcr, by = "ST_12")

# Read BCR map

bcr_shp = read_sf(paste0(database_dir, "bcr_strata.shp"))

# Read environmental data

Bioclim = read_csv(paste0(database_dir, "Bioclim.csv"))

# Trait labels, for plotting

trait_labels = c(`PC1_beak` = "Beak PC1", 
                 `PC1_wing` = "Wing PC1", 
                 `relative_bill_length` = "Relative Beak Length",
                 `relative_wing_length` = "Relative Wing Length", 
                 `Mass` = "Mass", 
                 `litter_or_clutch_size_n` = "Clutch Size", 
                 `GenLength` = "Generation Length", 
                 `log_Mass` = "log(Mass)", # after transformation
                 `log_clutch_size` = "log(Clutch Size)", # after transformation
                 `log_GenLength` = "log(Generation Length)", # after transformation
                 `corr_GenLength` = "Corrected Generation Length",
                 `relative_gen_length` = "Relative Generation Length"
)

color_4 = c("#2972B3", "#99C8F2", "#D98E04", "#734002")


########## Read PCA Results ##########


pca_loadings = read_csv(paste0(generated_dir, "PCA_Loadings.csv"))
pca_base = read_csv(paste0(generated_dir, "PCA_Base.csv"))
pca_topspp = read_csv(paste0(generated_dir, "PCA_TopSpp.csv"))
pca_restspp = read_csv(paste0(generated_dir, "PCA_RestSpp.csv"))

# Process data

pca_loadings = pca_loadings %>% 
  rename(Trait = `...1`)

pca_base_ind = pca_base %>% 
  rename("year_region" = "...1") %>% 
  separate_wider_delim(year_region, names = c("year", "region"), delim = "/", cols_remove = F) %>% 
  mutate(year = as.numeric(year))

pca_topspp_ind = pca_topspp %>% 
  rename("year_region" = "...1") %>% 
  separate_wider_delim(year_region, names = c("year", "region"), delim = "/", cols_remove = F) %>% 
  mutate(year = as.numeric(year))

pca_restspp_ind = pca_restspp %>% 
  rename("year_region" = "...1") %>% 
  separate_wider_delim(year_region, names = c("year", "region"), delim = "/", cols_remove = F) %>% 
  mutate(year = as.numeric(year))


########## Categorize Coordinates by Quadrants ##########


positions_base = pca_base_ind %>% 
  mutate(quad = case_when((PC1 > 0 & PC2 > 0) ~ 1,
                          (PC1 < 0 & PC2 > 0) ~ 2,
                          (PC1 < 0 & PC2 < 0) ~ 3,
                          (PC1 > 0 & PC2 < 0) ~ 4), 
         group = "base")

endpoints_base_start = positions_base %>% 
  filter(year == min(year, na.rm = T)) %>% 
  right_join(latlong_shp, by = c("region" = "ST_12")) %>% 
  mutate(year = replace_na(year, max(year, na.rm = T)))

endpoints_base_end = positions_base %>% 
  filter(year == max(year, na.rm = T)) %>% 
  right_join(latlong_shp, by = c("region" = "ST_12")) %>% 
  mutate(year = replace_na(year, min(year, na.rm = T)))

endpoints_base = rbind(endpoints_base_start, endpoints_base_end)

positions_topspp = pca_topspp_ind %>% 
  mutate(quad = case_when((PC1 > 0 & PC2 > 0) ~ 1,
                          (PC1 < 0 & PC2 > 0) ~ 2,
                          (PC1 < 0 & PC2 < 0) ~ 3,
                          (PC1 > 0 & PC2 < 0) ~ 4), 
         group = "topspp")

endpoints_topspp_start = positions_topspp %>% 
  filter(year == min(year, na.rm = T)) %>% 
  right_join(latlong_shp, by = c("region" = "ST_12")) %>% 
  mutate(year = replace_na(year, max(year, na.rm = T)))

endpoints_topspp_end = positions_topspp %>% 
  filter(year == max(year, na.rm = T)) %>% 
  right_join(latlong_shp, by = c("region" = "ST_12")) %>% 
  mutate(year = replace_na(year, min(year, na.rm = T)))

endpoints_topspp = rbind(endpoints_topspp_start, endpoints_topspp_end)

positions_restspp = pca_restspp_ind %>% 
  mutate(quad = case_when((PC1 > 0 & PC2 > 0) ~ 1,
                          (PC1 < 0 & PC2 > 0) ~ 2,
                          (PC1 < 0 & PC2 < 0) ~ 3,
                          (PC1 > 0 & PC2 < 0) ~ 4), 
         group = "restspp")

endpoints_restspp_start = positions_restspp %>% 
  filter(year == min(year, na.rm = T)) %>% 
  right_join(latlong_shp, by = c("region" = "ST_12")) %>% 
  mutate(year = replace_na(year, max(year, na.rm = T)))

endpoints_restspp_end = positions_restspp %>% 
  filter(year == max(year, na.rm = T)) %>% 
  right_join(latlong_shp, by = c("region" = "ST_12")) %>% 
  mutate(year = replace_na(year, min(year, na.rm = T)))

endpoints_restspp = rbind(endpoints_restspp_start, endpoints_restspp_end)

positions_all = rbind(positions_base, positions_topspp, positions_restspp)
endpoints_all = rbind(endpoints_base, endpoints_topspp, endpoints_restspp)

# Barplot for the overall study region

p_endpoints = 
  endpoints_all %>% 
  filter(!is.na(quad)) %>% 
  ggplot(aes(x = group, fill = as.factor(quad))) + 
  geom_bar(position = "fill", width = 0.5) + 
  scale_fill_manual(name = "Functional position (Quadrants)", values = color_4) +
  scale_x_discrete(name = "Focal Community", labels = c("base" = "Base", 
                                                        "topspp" = "Top", 
                                                        "restspp" = "Rest")) + 
  ylab("Number of Communities") + 
  facet_wrap(.~year) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        legend.position = "bottom")

p_endpoints

ggsave(paste0(figure_dir, "Endpoints.png"), p_endpoints, device = "png", 
       width = 1600, height = 1200, unit = "px")

# Functional positions over time in each region

tally_quad_time = positions_all %>% 
  select(year, region, quad, group) %>% 
  left_join(category_by_bcr, by = c("region" = "ST_12")) %>% 
  group_by(bioregion, group, year, quad) %>% 
  summarize(count = n()) %>% 
  pivot_wider(names_from = quad, values_from = count, values_fill = 0) %>% 
  rename(quad1 = `1`, quad2 = `2`, quad3 = `3`, quad4 = `4`) %>% 
  mutate(sum_grids = sum(quad1, quad2, quad3, quad4)) %>% 
  mutate(prop1 = quad1 / sum_grids, 
         prop2 = quad2 / sum_grids, 
         prop3 = quad3 / sum_grids, 
         prop4 = quad4 / sum_grids)


########## Visualize Trajectories ##########


# Plot single-grid example

example_grid_1 = pca_base_ind %>% 
  filter(region == "27_-98")

example_grid_2 = pca_base_ind %>% 
  filter(region == "40_-110")

scaling_factor = max(abs(example_grid_1$PC1), abs(example_grid_1$PC2), 
                     abs(example_grid_2$PC1), abs(example_grid_2$PC2)) / 
  max(abs(pca_loadings$PC1), abs(pca_loadings$PC2))
loadings_scaled = pca_loadings %>% 
  mutate(across(PC1:PC2, ~.*scaling_factor)) %>% 
  mutate(Trait = recode(Trait, !!!trait_labels))

p_pca = 
  ggplot() +
  geom_hline(yintercept = 0, color = "gray50") + 
  geom_vline(xintercept = 0, color = "gray50") +
  geom_segment(data = loadings_scaled, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm"))) +
  geom_point(data = example_grid_1, 
             aes(x = PC1, y = PC2, color = year)) + 
  geom_path(data = example_grid_1, 
             aes(x = PC1, y = PC2, color = year)) +
  geom_point(data = example_grid_2, 
             aes(x = PC1, y = PC2, color = year)) + 
  geom_path(data = example_grid_2, 
            aes(x = PC1, y = PC2, color = year)) +
  scale_color_gradientn(colors = color_4, name = "Year") + 
  xlim(-max(abs(loadings_scaled[, 2:3])), max(abs(loadings_scaled[, 2:3]))) +
  ylim(-max(abs(loadings_scaled[, 2:3])), max(abs(loadings_scaled[, 2:3]))) +
  xlab("PC1") + ylab("PC2") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

p_pca

ggsave(paste0(figure_dir, "Ind_Community_Vectors.png"), p_pca, device = "png", 
       width = 1600, height = 1200, unit = "px")

# Example grid, without PCA vectors

p_ind = 
  ggplot(example_grid_2, 
         aes(x = PC1, y = PC2, color = year)) + 
  geom_point() + 
  geom_path() + 
  geom_hline(yintercept = 0, color = "gray50") + 
  geom_vline(xintercept = 0, color = "gray50") +
  scale_color_gradientn(colors = color_4, name = "Year") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

p_ind

ggsave(paste0(figure_dir, "Ind_Community.png"), p_ind, device = "png", 
       width = 1600, height = 1200, unit = "px")


########## Calculate Trajectories and Metrics ##########


# Calculate trajectories

traj_base = calculate_traj(pca_base_ind)
traj_topspp = calculate_traj(pca_topspp_ind)
traj_restspp = calculate_traj(pca_restspp_ind)

# Calculate metrics of the temporal trajectories

traj_char_base = calculate_traj_metrics(traj_base)
traj_char_topspp = calculate_traj_metrics(traj_topspp)
traj_char_restspp = calculate_traj_metrics(traj_restspp)

# Extract values

mean_speed = tibble(base = traj_char_base$mean_speed, 
                    topspp = traj_char_topspp$mean_speed,
                    restspp = traj_char_restspp$mean_speed,
                    region = traj_char_base$region,
                    metric = "mean_speed")

Emax = tibble(base = traj_char_base$Emax, 
              topspp = traj_char_topspp$Emax, 
              restspp = traj_char_restspp$Emax,
              region = traj_char_base$region,
              metric = "Emax")

eff_displ = tibble(base = traj_char_base$eff_displ, 
                   topspp = traj_char_topspp$eff_displ, 
                   restspp = traj_char_restspp$eff_displ,
                   region = traj_char_base$region,
                   metric = "eff_displ")

eff_angle = tibble(base = traj_char_base$eff_angle, 
                   topspp = traj_char_topspp$eff_angle, 
                   restspp = traj_char_restspp$eff_angle,
                   region = traj_char_base$region,
                   metric = "eff_angle")

angle_bins = tibble(base = traj_char_base$angle_bins, 
                    topspp = traj_char_topspp$angle_bins,
                    restspp = traj_char_restspp$angle_bins,
                    region = traj_char_base$region,
                    metric = "angle_bins")

all_chars = rbind(mean_speed, Emax, eff_displ, eff_angle, angle_bins) %>% 
  pivot_longer(cols = base:restspp, names_to = "group", values_to = "value")


########## Plotting Maps of Grids ##########


# Combine communities with sf geometry

traj_char_base = traj_char_base %>% 
  right_join(latlong_shp, by = c("region" = "ST_12"))

traj_char_topspp = traj_char_topspp %>% 
  right_join(latlong_shp, by = c("region" = "ST_12"))

traj_char_restspp = traj_char_restspp %>% 
  right_join(latlong_shp, by = c("region" = "ST_12"))

# Range of values to unify gradients across maps
# Mean speed and effective displacement are sqrt-transformed

range_mean_speed = sqrt(range(all_chars %>% filter(metric == "mean_speed") %>% pull(value), na.rm = T))
range_Emax = range(all_chars %>% filter(metric == "Emax") %>% pull(value), na.rm = T)
range_eff_displ = sqrt(range(all_chars %>% filter(metric == "eff_displ") %>% pull(value), na.rm = T))

# Maps, base communities

pls_base = plot_maps(endpoints_base, traj_char_base, comm_type = "base")
pls_base

# Maps, top contributors

pls_topspp = plot_maps(endpoints_topspp, traj_char_topspp, comm_type = "topspp")
pls_topspp

# Maps, rest of the communities

pls_restspp = plot_maps(endpoints_restspp, traj_char_restspp, comm_type = "restspp")
pls_restspp

p_start_end_angle_all = pls_base$start_end_angle + pls_topspp$start_end_angle + pls_restspp$start_end_angle

ggsave(paste0(figure_dir, "Start_End_Angle_all.png"), p_start_end_angle_all, device = "png", 
       width = 4000, height = 3000, unit = "px")


########## Spatial Autocorrelation ##########


# Calculate spatial heterogeneity by Moran's I

sp_hetero_base = calculate_spatial_heterogeneity(positions_base, latlong_shp, "base")
sp_hetero_topspp = calculate_spatial_heterogeneity(positions_topspp, latlong_shp, "topspp")
sp_hetero_restspp = calculate_spatial_heterogeneity(positions_restspp, latlong_shp, "restspp")

# Stats

sp_hetero_all = rbind(sp_hetero_base, sp_hetero_topspp, sp_hetero_restspp)

summary(aov(`Moran I statistic` ~ year + community, data = sp_hetero_all))

broom::tidy(lm(`Moran I statistic` ~ year, data = sp_hetero_base))
broom::tidy(lm(`Moran I statistic` ~ year, data = sp_hetero_topspp))
broom::tidy(lm(`Moran I statistic` ~ year, data = sp_hetero_restspp))

# Visualization

p_hetero = sp_hetero_all %>% 
  ggplot(aes(x = year, y = `Moran I statistic`, color = community, fill = community, shape = community)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_hline(yintercept = 0) + 
  scale_color_manual(name = "Community", labels = c("Base", "Rest", "Top"), 
                     values = color_4[1:3]) +
  scale_fill_manual(name = "Community", labels = c("Base", "Rest", "Top"), 
                     values = color_4[1:3]) + 
  scale_shape_manual(name = "Community", labels = c("Base", "Rest", "Top"), 
                     values = c(16, 17, 15)) + 
  xlab("Year") + ylab("Spatial heterogeneity (Moran's I)") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10), 
        legend.position = "bottom")

p_hetero

ggsave(paste0(figure_dir, "Spatial_Heterogeneity.png"), p_hetero, device = "png", 
       width = 1600, height = 1200, unit = "px")


########## Tally Values by BCR ##########


# Endpoints and effective angle barplots by bioregions

end_angle_region = endpoints_all %>% 
  filter(!is.na(quad)) %>% 
  select(region, group, quad, year) %>% 
  mutate(year = as.character(year)) %>% 
  bind_rows(all_chars %>% 
         filter(metric == "angle_bins") %>% 
         mutate(year = "Effective Angle") %>% 
         rename(quad = value) %>% 
         select(-metric)) %>% 
  left_join(category_by_bcr, by = c("region" = "ST_12"))

p_end_angle_bio = end_angle_region %>% 
  ggplot(aes(x = year, fill = as.factor(quad))) + 
  geom_bar(position = "fill", width = 0.5) + 
  scale_fill_manual(name = "Functional Position (Quadrants)", values = color_4) +
  xlab("Start (1970), End (2021), and Effective Angle") + 
  ylab("Proportion of Communities") + 
  facet_grid(bioregion~group, labeller = labeller(group = c("base" = "Base", 
                                                            "topspp" = "Top", 
                                                            "restspp" = "Rest"))) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        legend.position = "bottom")

ggsave(paste0(figure_dir, "Endpoints_Bio_Angles.png"), p_end_angle_bio, device = "png", 
       width = 2400, height = 3200, unit = "px")

# BCR maps by bioregions

p_bcr = bcr_shp %>% 
  left_join(category_by_bcr, by = "strt_nm") %>% 
  ggplot() + 
  geom_sf(aes(fill = bioregion, geometry = geometry)) +
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        legend.position = "bottom", 
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 10))

ggsave(paste0(figure_dir, "Bioregions.png"), p_bcr, device = "png", 
       width = 2400, height = 1800, unit = "px")


########## Quadrant-Environment Associations ##########


# Compare temperature seasonality by effective angle bins, 2 and 3 only

angle_bins_env = all_chars %>% 
  filter(metric == "angle_bins") %>% 
  filter(value == 2 | value == 3) %>% 
  left_join(Bioclim, by = c("region" = "ST_12"), relationship = "many-to-many") %>% 
  mutate(temp_seasonality = bio5 - bio6) %>% 
  group_by(region, value, group) %>% 
  summarize(mean_temp_seasonality = mean(temp_seasonality))

p_angles_env = angle_bins_env %>% 
  ggplot(aes(x = as.factor(value), y = mean_temp_seasonality, 
             fill = as.factor(value))) + 
  geom_boxplot() + 
  ggpubr::stat_compare_means(method = "t.test") + 
  xlab("Effective Angle") + ylab("Mean Temperature Seasonality over Years") +
  scale_fill_manual(name = "Effective Angle", values = c("#99C8F2", "#D98E04")) +
  facet_wrap(.~group, labeller = labeller(group = c("base" = "Base", 
                                                    "topspp" = "Top", 
                                                    "restspp" = "Rest"))) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12), 
        legend.position = "none")

ggsave(paste0(figure_dir, "Angles_TempSeasonality.png"), p_angles_env, device = "png", 
       width = 1600, height = 1200, unit = "px")
