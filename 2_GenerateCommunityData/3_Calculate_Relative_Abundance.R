# CALCULATE RELATIVE ABUNDANCE OF EACH SPECIES IN EACH GRID
# Hengxing Zou


########## Initialization ##########


library(tidyverse)

# Change the directory according to your working environment

generated_dir = "GeneratedData/"


########## Calculation ##########


# Read data

All_Comms = read_csv(paste0(generated_dir, "All_Communities_Filtered.csv"))

CWM_All = read_csv(paste0("CWMetrics_Filtered.csv"))

# Calculate relative abundance of species in each grid and year for the whole community

comm_abundance = CWM_All %>% 
  select(Region, year, Community_Abundance) %>% 
  distinct()

relative_abundance = All_Comms %>% 
  left_join(comm_abundance, by = c("region" = "Region", "year")) %>% 
  mutate(Relative_Abundance = index / Community_Abundance)

# Write relative abundance data

write_csv(relative_abundance, paste0(generated_dir, "Relative_Abundance.csv"))
