# CALCULATE SPECIES CONTRIBUTIONS TO TEMPORAL TRENDS OF CWM
# Hengxing Zou


########## Initialization ##########


library(tidyverse)
library(foreach)
library(doSNOW)

cl = makeCluster(4, outfile = "")
registerDoSNOW(cl)

# Change the directory according to your working environment

generated_dir = "GeneratedData/"
database_dir = "Database/"


########## Read and Process Data ##########


# Read community-weighted means of all species

All_Weighted_Metrics = read_csv(paste0(generated_dir, "CWMetrics_Filtered.csv")) %>% 
  rename("Hand.Wing.Index" = "Hand-Wing.Index")

# Join CWM with climate and anthrome data

cwm_mean = All_Weighted_Metrics %>% 
  filter(Metric == "CWM") %>% 
  mutate(year_since_1969 = year - 1969)

# Read relative abundance data

Relative_Abundance = read_csv(paste0(generated_dir, "Relative_Abundance_orig.csv"))

rel_abundance = Relative_Abundance %>% 
  select(region, latin, year, index, Community_Abundance, Relative_Abundance) %>% 
  mutate(year_since_1969 = year - 1969)

# Read functional trait data

All_Funct_Traits = read_csv(paste0(database_dir, "Filtered_Funct_Data.csv"))

selected_traits = All_Funct_Traits %>% 
  select(Species, PC1_beak, PC1_wing, relative_wing_length, relative_bill_length, litter_or_clutch_size_n, corr_GenLength, Mass) %>% 
  column_to_rownames("Species")

md_traits = colnames(selected_traits)

# Remove big original datasets

rm(All_Weighted_Metrics, Relative_Abundance, All_Funct_Traits)


########## Calculate Contributions of Each Species ##########


grids = unique(cwm_mean$Region)

# Compute slopes over time

all_contribs = foreach(g = grids, .combine = rbind, .packages = "tidyverse") %dopar% {
  
  cat("Calculating temporal slopes for", g, "\n")
  
  cwm_mean_grid = cwm_mean %>% 
    filter(Region == g)
  rel_abun_grid = rel_abundance %>% 
    filter(region == g)
  spp_list = rel_abun_grid %>% 
    pull(latin) %>% 
    unique()
    
  temporal_var = var(cwm_mean_grid$year_since_1969)
  
  # Calculate slopes of relative abundances for each species
  slopes_spp = numeric()
  
  for (spp in spp_list) {
    
    rel_abun_spp = rel_abun_grid %>% 
      filter(latin == spp) %>% 
      select(year_since_1969, Relative_Abundance)
    
    # For species with data-deficient years, add 0 to missing years
    if (length(rel_abun_spp$year_since_1969) < length(cwm_mean_grid$year_since_1969)) {
      
      missing_years = setdiff(cwm_mean_grid$year_since_1969, rel_abun_spp$year_since_1969)
      rel_abun_spp = rbind(rel_abun_spp, tibble(year_since_1969 = missing_years, 
                                                Relative_Abundance = 0)) %>% 
        arrange(year_since_1969)
      
    }
    
    slopes_spp[spp] = cov(cwm_mean_grid$year_since_1969, rel_abun_spp$Relative_Abundance) / temporal_var
    
  }
  
  # Calculate contributions of each species to CWM
  slopes_cwm = numeric()
  output_tr = tibble(.rows = length(spp_list))
  
  for (tr in md_traits) {
    
    # Calculate linear slope of each CWM
    slope_cwm = cov(cwm_mean_grid$year_since_1969, cwm_mean_grid[, tr]) / temporal_var
    slopes_cwm = c(slopes_cwm, slope_cwm)
    
    # Calculate the contribution of each species to the overall CWM slope
    contrib_spp = slopes_spp * selected_traits[spp_list, tr]
    output_tr = cbind(output_tr, contrib_spp)
  
  }
  
  # Record species abundance slopes
  abun_slopes = as.data.frame(slopes_spp) %>% 
    rownames_to_column("Species") %>% 
    rename(Abun_Slope = slopes_spp)

  # Record species contributions
  colnames(output_tr) = md_traits
  contributions = output_tr %>% 
    rownames_to_column("Species") %>% 
    pivot_longer(-Species, names_to = "Trait", values_to = "Species_Contrib")
  
  # Record CWM ~ Time slopes and combine with all others
  names(slopes_cwm) = md_traits
  output = as.data.frame((slopes_cwm)) %>% 
    rownames_to_column("Trait") %>% 
    rename(CWM_Slope = `(slopes_cwm)`) %>% 
    full_join(contributions, by = "Trait") %>% 
    full_join(abun_slopes, by = "Species") %>% 
    mutate(Region = g)
  
  output
  
}

# Write data

write_csv(all_contribs, paste0(generated_dir, "Species_Contributions.csv"))
