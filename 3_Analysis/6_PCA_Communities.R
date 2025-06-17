# CONSTRUCTING THE FUNCTIONAL SPACE FOR ALL COMMUNITIES
# Hengxing Zou


########## Initialization ##########


# Double check to make sure the directory is correct

source("3_Analysis/1_Read_All_Data.R")

library(factoextra)

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


########## PCA of Base Community ##########


# Convert each community (region-year) into a single row

cwm_mean_rownames = cwm_mean %>% 
  unite(year_region, c("year", "Region"), sep = "/") %>% 
  select(1, 12, 14, 16:19, 21) %>% 
  column_to_rownames(var = "year_region")

# Do PCA for the full community

pca_base = prcomp(cwm_mean_rownames, scale = T, center = T)
summary(pca_base)

# Check axes (quick visualization)

fviz_pca_var(pca_base)

# Flip axes such that relative beak/ wing length arrows are positive for PC1, 
# and clutch size arrow is positive for PC2

pca_base$rotation[, 1:2] = -pca_base$rotation[, 1:2]
pca_base$x[, 1:2] = -pca_base$x[, 1:2]

# Extract axes for further analyses

loadings_base = pca_base$rotation

# Extract PCA scores of communities

coord_base = pca_base$x

# Percentage variance explained

variances = (pca_base$sdev)^2
prop_explained = variances / sum(variances) * 100

# Better visualization

p_pca = as.data.frame(loadings_base) %>%
  rownames_to_column("Trait") %>% 
  mutate(Trait = recode(Trait, !!!trait_labels)) %>% 
  ggplot() +
  geom_hline(yintercept = 0, color = "gray50") + 
  geom_vline(xintercept = 0, color = "gray50") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm"))) +
  ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = Trait),
                           hjust = ifelse(loadings_base[, 1] > 0, -0.1, 1.2), vjust = 1.4, 
                           segment.color = "transparent", 
                           size = 3) + 
  annotate(geom = "text", x = c(0.5, -0.5, -0.5, 0.5), y = c(0.5, 0.5, -0.5, -0.5), 
           label = c("Quadrant 1", "Quadrant 2", "Quadrant 3", "Quadrant 4"), 
           size = 6, color = color_4) +
  xlim(-max(abs(loadings_base[, 1:2])), max(abs(loadings_base[, 1:2]))) + 
  ylim(-max(abs(loadings_base[, 1:2])), max(abs(loadings_base[, 1:2]))) + 
  xlab(paste0("PC1 (", round(prop_explained[1], 2), "%)")) + 
  ylab(paste0("PC2 (", round(prop_explained[2], 2), "%)")) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12))

p_pca

# Save figure

ggsave(paste0(figure_dir, "Functional_Space_Axes.png"), p_pca, device = "png", 
       width = 1600, height = 1200, unit = "px")


########## PCA of Communities by Trends ##########


# Extract communities by trends

cwm_mean_topspp_rownames = cwm_mean_topspp %>% 
  unite(year_region, c("year", "Region"), sep = "/") %>% 
  select(1, 12, 14, 16:20) %>% 
  column_to_rownames(var = "year_region")

cwm_mean_restspp_rownames = cwm_mean_restspp %>% 
  unite(year_region, c("year", "Region"), sep = "/") %>% 
  select(1, 12, 14, 16:20) %>% 
  column_to_rownames(var = "year_region")

# Center data based on the base PCA

cwm_mean_topspp_rownames_centered = sweep(cwm_mean_topspp_rownames, 2, pca_base$center)
cwm_mean_topspp_rownames_scales = as.matrix(sweep(cwm_mean_topspp_rownames_centered, 2, pca_base$scale, "/"))

cwm_mean_restspp_rownames_centered = sweep(cwm_mean_restspp_rownames, 2, pca_base$center)
cwm_mean_restspp_rownames_scales = as.matrix(sweep(cwm_mean_restspp_rownames_centered, 2, pca_base$scale, "/"))

# Project data onto the base principle components

pca_topspp_scores = cwm_mean_topspp_rownames_scales %*% loadings_base
pca_restspp_scores = cwm_mean_restspp_rownames_scales %*% loadings_base


########### Write Data ##########


write.csv(loadings_base, paste0(generated_dir, "PCA_Loadings.csv"))
write.csv(coord_base, paste0(generated_dir, "PCA_Base.csv"))
write.csv(pca_topspp_scores, paste0(generated_dir, "PCA_TopSpp.csv"))
write.csv(pca_restspp_scores, paste0(generated_dir, "PCA_RestSpp.csv"))

