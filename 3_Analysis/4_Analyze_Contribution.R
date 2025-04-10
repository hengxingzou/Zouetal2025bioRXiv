# ANALYZE SPECIES CONTRIBUTIONS TO TEMPORAL TRENDS IN CWM
# Hengxing Zou


########## Initialization ##########


library(tidyverse)
library(patchwork)
library(sf)

library(foreach)
library(doSNOW)

cl = makeCluster(4, outfile = "")
registerDoSNOW(cl)

color_2 = c("#2972B3", "#D98E04")

# Trait labels for plotting

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
                 `corr_GenLength` = "Corrected Generation Length"
)

# Change the directory according to your working environment

generated_dir = "GeneratedData/"
database_dir = "Database/"
figure_dir = "Figures/"


########## Read Data ##########


# Read functional trait data

All_Funct_Traits = read_csv(paste0(database_dir, "../Databases/Filtered_Funct_Data.csv"))

selected_traits = All_Funct_Traits %>% 
  select(Species, PC1_beak, PC1_wing, relative_wing_length, relative_bill_length, litter_or_clutch_size_n, corr_GenLength, Mass) %>% 
  column_to_rownames("Species")

md_traits = colnames(selected_traits)

selected_tr_scaled = selected_traits %>% 
  mutate(across(!c(relative_wing_length, relative_bill_length), ~scale(., scale = T, center = T))) %>% 
  mutate_all(~str_remove(., "[, 1]")) %>%
  mutate_all(~as.numeric(.)) %>% 
  rownames_to_column("Species") %>% 
  pivot_longer(2:8, names_to = "Trait", values_to = "Value")

# Read species contributions

All_Contribs = read_csv(paste0(generated_dir, "Species_Contributions.csv"))

# Read map

latlong_shp = read_sf(paste0(database_dir, "BBS_LatLong_strata.shp"))


########## Top 10% Contributing Species ##########


# All species, with scaled functional traits

all_contribs_tr = All_Contribs %>% 
  select(-CWM_Slope) %>% 
  
  # calculate median slope of abundance by species
  group_by(Species) %>% 
  mutate(median_slope = median(Abun_Slope), 
         q5 = quantile(Abun_Slope, 0.05), 
         q95 = quantile(Abun_Slope, 0.95)) %>% 
  
  group_by(Region, Trait) %>% 
  left_join(selected_tr_scaled, by = c("Species", "Trait")) %>% 
  mutate(Contrib_Sign = if_else(Species_Contrib >= 0, "pos", "neg")) %>% 
  rename(Scaled_Trait_Val = Value)

# Rank species by contribution, then filter by top 10%

top_contribs = all_contribs_tr %>% 
  group_by(Region, Trait) %>% 
  mutate(threshold = quantile(abs(Species_Contrib), 0.9)) %>% 
  filter(abs(Species_Contrib) >= threshold) %>% 
  select(-threshold)

# Tally number of grids in which the species is a top contributor

count_grids = top_contribs %>% 
  group_by(Trait, Species, Contrib_Sign) %>% 
  summarize(Num_Top_Grids = n_distinct(Region)) %>% 
  group_by(Trait, Species) %>% 
  mutate(Num_Total_Top_Grids = sum(Num_Top_Grids)) %>% 
  group_by(Trait)

# Find the top 10 species

top_10_signs = count_grids %>% 
  
  # get the top 10 for further analysis
  # need to use n = 20 because Num_Total_Top_Grids has two identical values for each species
  # this is more conservative (less species)
  slice_max(order_by = Num_Total_Top_Grids, n = 20, with_ties = T)
  
  # get the top 10%
  # mutate(threshold = quantile(Num_Total_Top_Grids, 0.9)) %>%
  # filter(Num_Total_Top_Grids >= threshold) %>%
  # select(-threshold)
  
# This species list is used to calculate CWMs as the top contributors

top_10_spp = unique(top_10_signs$Species)

# This data frame contains more info for each species, e.g., median abundance slope

top_10_tr = top_10_signs %>% 
  select(-Contrib_Sign, -Num_Top_Grids) %>% 
  distinct() %>% 
  left_join(top_contribs %>% 
              ungroup() %>% 
              select(-Species_Contrib, -Abun_Slope, -Region, -Contrib_Sign) %>% 
              distinct(), 
            by = c("Species", "Trait"))

# Optional: write data for further analysis

# write_csv(top_10_tr, paste0(generated_dir, "Top_10_Traits.csv"))

# Visualization

pls_count = list()

for (i in 1:length(md_traits)) {
  
  pl_ct = top_10_signs %>% 
    filter(Trait == md_traits[i]) %>% 
    ggplot(aes(x = Num_Top_Grids, y = reorder(Species, Num_Total_Top_Grids), fill = Contrib_Sign)) + 
    ggtitle(recode(md_traits[i], !!!trait_labels)) +
    geom_bar(stat = "identity") + 
    scale_fill_manual(name = "Sign of contribution", values = color_2, labels = c("Negative", "Positive")) + 
    xlab("Number of Grids") + ylab("Dominant Species") + 
    theme_bw() + 
    theme(title = element_text(size = 12), 
          axis.text.y = element_text(size = 10, angle = 30, hjust = 0.75),
          axis.text.x = element_text(size = 10),
          axis.title = element_text(size = 12), 
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.position = "bottom"
    )
  
  pls_count[[i]] = pl_ct
  
}

names(pls_count) = md_traits

pls_count

# Save figures

ggsave(paste0(figure_dir, "Contributions_Beak.png"), pls_count$PC1_beak, device = "png", 
       width = 1600, height = 1200, unit = "px")
ggsave(paste0(figure_dir, "Contributions_Wing.png"), pls_count$PC1_wing, device = "png", 
       width = 1600, height = 1200, unit = "px")
ggsave(paste0(figure_dir, "Contributions_Relative_Wing_Length.png"), pls_count$relative_wing_length, device = "png", 
       width = 1600, height = 1200, unit = "px")
ggsave(paste0(figure_dir, "Contributions_Relative_Beak_Length.png"), pls_count$relative_bill_length, device = "png", 
       width = 1600, height = 1200, unit = "px")
ggsave(paste0(figure_dir, "Contributions_Clutch_Size.png"), pls_count$litter_or_clutch_size_n, device = "png", 
       width = 1600, height = 1200, unit = "px")
ggsave(paste0(figure_dir, "Contributions_Corr_Generation_Length.png"), pls_count$corr_GenLength, device = "png", 
       width = 1600, height = 1200, unit = "px")
ggsave(paste0(figure_dir, "Contributions_Mass.png"), pls_count$Mass, device = "png", 
       width = 1600, height = 1200, unit = "px")


########## Analyze Top Contributors, Contour Plots ##########


# Optional: read top contributors from saved data

# top_10_tr = read_csv(paste0(generated_dir, "Top_10_Traits.csv"))

# Visualization

pls_avg_topspp = list()

for (i in 1:length(md_traits)) {

  top_10_spp = top_10_tr %>% 
    filter(Trait == md_traits[i]) %>% 
    pull(Species)
  
  all_spps = count_grids %>% 
    filter(Trait == md_traits[i]) %>% 
    select(-Num_Top_Grids, -Contrib_Sign) %>% 
    distinct() %>% 
    left_join(top_contribs %>% 
                ungroup() %>% 
                select(-Species_Contrib, -Abun_Slope, -Region, -Contrib_Sign) %>% 
                distinct(), 
              by = c("Species", "Trait")) %>% 
    mutate(Top_Contrib = if_else(Species %in% top_10_spp, "yes", "no"))

  xrange = range(all_spps$Scaled_Trait_Val)*1.05
  yrange = range(all_spps$median_slope)*1.05
  # yrange = range(all_spps$q5, tr_dat$q95)*1.05
  cont = expand_grid(x = seq(xrange[1], xrange[2], by = (xrange[2]-xrange[1])/150), 
                     y = seq(yrange[1], yrange[2], by = (yrange[2]-yrange[1])/150)) %>% 
    mutate(z = x*y)
  zrange = range(cont$z)
  
  pl_avg_topspp = ggplot() +
    geom_tile(data = cont, aes(x = x, y = y, fill = z)) +
    geom_contour(data = cont, aes(x = x, y = y, z = z), color = "gray50", 
                 breaks = setdiff(seq(zrange[1], zrange[2], length = 10), 0)) + 
    geom_point(data = all_spps,
               aes(x = Scaled_Trait_Val, y = median_slope,
                   alpha = Top_Contrib, size = Num_Total_Top_Grids)) +
    # geom_pointrange(data = all_spps,
    #                 aes(x = Scaled_Trait_Val, y = median_slope, ymin = q5, ymax = q95,
    #                     alpha = Top_Contrib)) +
    ggrepel::geom_text_repel(data = all_spps %>% filter(Top_Contrib == "yes"), 
                             aes(x = Scaled_Trait_Val, y = median_slope, label = Species)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray25") + 
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray25") +
    coord_cartesian() +
    ggtitle(recode(md_traits[i], !!!trait_labels)) +
    xlab("Scaled Trait Value") + 
    ylab("Median of Abundance Slopes") + 
    scale_fill_gradient2(low = color_2[1], mid = "white", high = color_2[2], midpoint = 0,
                         name = "Species Contribution", 
                         breaks = c(zrange[1], 0, zrange[2]), 
                         labels = c("-", "0", "+")) +
    scale_size(name = "Number of Grids where \nSpecies rank Top 10% \nby Contributions") +
    scale_alpha_manual(values = c(0.1, 1), guide = "none") +
    theme_bw() + 
    theme(title = element_text(size = 12), 
          panel.grid = element_blank(), 
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10), 
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.position = "bottom")
  
  pls_avg_topspp[[i]] = pl_avg_topspp
  
}

names(pls_avg_topspp) = md_traits

pls_avg_topspp

# Save single-panel figures

ggsave(paste0(figure_dir, "Beak_Contour.png"), pls_avg_topspp$PC1_beak, 
       width = 3200, height = 2400, unit = "px")
ggsave(paste0(figure_dir, "Wing_Contour.png"), pls_avg_topspp$PC1_wing, 
       width = 3200, height = 2400, unit = "px")
ggsave(paste0(figure_dir, "Relative_Wing_Length_Contour.png"), pls_avg_topspp$relative_wing_length, 
       width = 3200, height = 2400, unit = "px")
ggsave(paste0(figure_dir, "Relative_Bill_Length_Contour.png"), pls_avg_topspp$relative_bill_length, 
       width = 3200, height = 2400, unit = "px")
ggsave(paste0(figure_dir, "Clutch_Size_Contour.png"), pls_avg_topspp$litter_or_clutch_size_n, 
       width = 3200, height = 2400, unit = "px")
ggsave(paste0(figure_dir, "Corr_GenLength_Contour.png"), pls_avg_topspp$corr_GenLength, 
       width = 3200, height = 2400, unit = "px")
ggsave(paste0(figure_dir, "Mass_Contour.png"), pls_avg_topspp$Mass, 
       width = 3200, height = 2400, unit = "px")

# Combine all panels

pls_contour = (pls_avg_topspp$PC1_beak + theme(legend.position = "none")) + 
  (pls_avg_topspp$PC1_wing + theme(legend.position = "none")) + 
  (pls_avg_topspp$relative_bill_length + theme(legend.position = "none")) + 
  (pls_avg_topspp$relative_wing_length + theme(legend.position = "none")) + 
  (pls_avg_topspp$Mass + theme(legend.position = "none")) + 
  (pls_avg_topspp$litter_or_clutch_size_n + theme(legend.position = "none")) + 
  (pls_avg_topspp$corr_GenLength + theme(legend.position = "none")) + 
  plot_layout(ncol = 4, nrow = 2, axes = "collect")

ggsave(paste0(figure_dir, "All_Contours.png"), pls_contour, 
       width = 6400, height = 3200, unit = "px")
