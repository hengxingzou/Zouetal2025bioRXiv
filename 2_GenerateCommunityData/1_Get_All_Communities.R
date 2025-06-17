# COMBINE ALL SPECIES INDICES INTO COMMUNITIES
# Hengxing Zou


########## Initialization ##########


library(tidyverse)
library(bbsBayes2)
library(foreach)
library(doSNOW)

cl = makeCluster(8, outfile = "")
registerDoSNOW(cl)

# Change the directory according to your working environment

indices_dir = "Indices/"
database_dir = "Database/"
generated_dir = "GeneratedData/"


########## Read All Indices ##########


# Read all indices from model fitting

files_indices = (Sys.glob(paste0(indices_dir, "*.csv")))

all_indices_lst = lapply(files_indices, function(x) read_csv(x))

All_Indices = bind_rows(all_indices_lst)

# Optional: write all indices

# write_csv(All_Indices, paste0(indices_dir, "All_Indices.csv"))


########## Generate Abundances for All Communities ##########


# Optional: read saved combined indices

# All_Indices = read_csv(paste0(indices_dir, "All_Indices.csv"))

# Read species list, annotated by number of records, guild, and model fitting diagnostics

All_Spp_List = read_csv(paste0(database_dir, "BBS_List_PIF"))

# Remove duplicates of common names due to mergers

spp_list_unique = All_Spp_List %>% 
  select(-aou, -num_records) %>% 
  distinct()

# Find all grids

grids = All_Indices %>% 
  filter(region_type != "continent") %>% 
  pull(region) %>% 
  unique()

# Get all communities

All_Communities = foreach(g = grids, .combine = rbind, .packages = "tidyverse") %dopar% {
  
  cat("Calculating community matrix for the grid", g, "\n")
  
  local_comm = All_Indices %>%
    filter(region == g)
  
  comm_mat = local_comm %>% 
    select(year, index, index_q_0.05, index_q_0.95, species) %>% 
    rename(english = species) %>% 
    left_join(spp_list_unique, by = "english") %>% 
    select(-english, -`exclude?`, -`rerun`, -`convergence`) %>%
    mutate(region = g)
  
  comm_mat
  
}

# Optional: write all unfiltered communities

# write_csv(All_Communities, paste0(generated_dir, "All_Communities.csv"))


########## Filter Communities and Species ##########


# Communities (grids) and species are filtered after fitting all the models in case these filters need to be changed
# in further analyses

# Optional: read unfiltered communities from saved file

# All_Communities = read_csv(generated_dir, "All_Communities.csv")

# Read functional traits

All_Funct_Data = read_csv(paste0(database_dir, "All_Funct_Data.csv")) %>% 
  column_to_rownames(var = "Species")

# Remove nocturnal and pelagic species according to EltonTraits

Nocturnal_Pelagic = All_Funct_Data %>% 
  select(PelagicSpecialist, Nocturnal) %>% 
  filter(PelagicSpecialist == 1 | Nocturnal == 1)

All_Comms_Filtered = All_Comms %>% 
  filter(!latin %in% rownames(Nocturnal_Pelagic)) %>% 
  # same year range as trends data
  filter(year >= 1970, year <= 2021) %>% 
  separate_wider_delim(region, names = c("lat", "long"), delim = "_", cols_remove = F) %>% 
  mutate(across(c("lat", "long"), as.numeric)) %>% 
  filter(lat <= 55, long >= -135)

# Filter communities with less than 10 species in any year

comms_remove = All_Comms_Filtered %>% 
  group_by(region, year) %>% 
  summarize(nspp = n_distinct(latin)) %>% 
  filter(nspp <= 10)

All_Comms_Filtered = All_Comms_Filtered %>% 
  filter(!(region %in% comms_remove$region))

# Apply PIF adjustments

All_Comms_Filtered_adj = All_Comms_Filtered %>% 
  left_join(All_Spp_List %>% select(-aou, -english) %>% distinct(), by = "latin") %>% 
  mutate(index = index * PIF_adj_birds_km2, 
         index_q_0.05 = index_q_0.05 * PIF_adj_birds_km2, 
         index_q_0.95 = index_q_0.95 * PIF_adj_birds_km2) %>% 
  select(-PIF_adj_birds_km2)


########## Write Data ##########


# Write the filtered species list

filter_spp_lst = All_Spp_List %>% 
  filter(latin %in% All_Comms_Filtered$latin)

write_csv(filter_spp_lst, paste0(database_dir, "Filtered_Spp_List.csv"))

# Write the filtered community data for future use

write_csv(All_Comms_Filtered, paste0(generated_dir, "All_Communities_Filtered.csv"))
