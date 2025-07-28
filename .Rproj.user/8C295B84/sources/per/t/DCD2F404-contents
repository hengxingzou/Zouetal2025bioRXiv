# CALCULATE COMMUNITY-WEIGHTED MEANS FOR ALL OR SELECTED SPECIES
# Hengxing Zou


########## Initialization ##########


library(tidyverse)
library(foreach)
library(doSNOW)

cl = makeCluster(10, outfile = "")
registerDoSNOW(cl)

# Double check to make sure the directory is correct

source("2_GenerateCommunityData/calculate_funct_distr.R")

# Change the directory according to your working environment

database_dir = "Database/"
generated_dir = "GeneratedData/"


########## Read Data ##########


# Read all communities with trends; 
# All_Communities_Filtered is the one that filters out communities where 
# species with either positive or negative trends are less than 10

All_Comms = read_csv(paste0(database_dir, "All_Communities_Filtered.csv"))

top_spp = c("Branta canadensis", "Turdus migratorius", "Cathartes aura", "Zenaida macroura",
            "Meleagris gallopavo", "Corvus brachyrhynchos", "Quiscalus quiscula", "Corvus corax",
            "Sturnella neglecta", "Phasianus colchicus", "Passer domesticus", "Spizella passerina",
            "Petrochelidon pyrrhonota", "Polioptila caerulea", "Vireo olivaceus", "Eremophila alpestris",
            "Passerculus sandwichensis", "Melospiza melodia", "Cardinalis cardinalis", "Troglodytes aedon",
            "Geothlypis trichas", "Sturnus vulgaris", "Agelaius phoeniceus", "Molothrus ater", "Colinus virginianus")

# IMPORTANT: TO CALCULATE CWM OF ONLY TOP SPECIES, UNCOMMENT BELOW

# All_Comms = All_Comms %>%
#   filter(latin %in% top_spp)

# IMPORTANT: TO CALCULATE CWM OF THE REST OF THE SPECIES, UNCOMMENT BELOW

# All_Comms = All_Comms %>%
#   filter(!(latin %in% top_spp))

# Read functional traits

All_Funct_Data = read_csv(paste0(database_dir, "Filtered_Funct_Data.csv")) %>% 
  column_to_rownames(var = "Species")

# Select functional traits

traits_exten = All_Funct_Data[, c(1:10, 60:65, # morphology
                                  11, 38, 59, 66:67)] # life history

grids = unique(All_Comms$region)


########## Calculate CWM for All Species ##########


weighted = foreach (g = grids, .combine = rbind, .packages = "tidyverse") %dopar% {
  
  cat("Calculating community-weighted metrics for", g, "\n")
                     
  local_comm = All_Comms %>% 
    filter(region == g) %>% 
    select(-starts_with("trend"), -percent_change)
  
  comm_mat = local_comm %>% 
    select(-index_q_0.05, -index_q_0.95, -region, -lat, -long) %>% 
    pivot_wider(names_from = latin, values_from = index) %>% 
    replace(is.na(.), 0) %>%
    column_to_rownames(var = "year")
  
  cw_metrics = calculate_funct_distr(as.matrix(traits_exten[colnames(comm_mat), ]), 
                                     as.matrix(comm_mat))
  
  cwm = cw_metrics[[1]] %>% 
    mutate(Metric = "CWM", Region = g, Index = "Mean", 
           Community_Abundance = rowSums(comm_mat)) %>% 
    rownames_to_column(var = "year")
  
  cwv = cw_metrics[[2]] %>% 
    mutate(Metric = "CWV", Region = g, Index = "Variance", 
           Community_Abundance = rowSums(comm_mat)) %>% 
    rownames_to_column(var = "year")
  
  grid_results = rbind(cwm, cwv) %>% 
    mutate(year = as.numeric(year)) %>% 
    arrange(year) %>% 
    separate_wider_delim(Region, names = c("lat", "long"), delim = "_", cols_remove = F) %>% 
    mutate(across(c("lat", "long"), as.numeric))
  
  grid_results
  
}

# Write data

write_csv(weighted, paste0(generated_dir, "CWMetrics_Filtered.csv"))

# IMPORTANT: TO WRITE DATA FOR TOP SPECIES ONLY, UNCOMMENT BELOW

write_csv(weighted, paste0(generated_dir, "CWMetrics_TopSpp.csv"))

# IMPORTANT: TO WRITE DATA FOR THE REST OF THE SPECIES, UNCOMMENT BELOW

write_csv(weighted, paste0(generated_dir, "CWMetrics_RestSpp.csv"))
