# RUN MODELS TO ESTIMATE SPECIES ABUNDANCE
# Hengxing Zou

# WARNING: RUNNING MODELS CAN TAKE A VERY LONG TIME, EVEN DAYS
# THIS VERSION OF THE SCRIPT IS FOR RUNNING ON HIGH-PERFORMANCE CLUSTERS
# ADJUST CODE AND DIRECTORIES ACCORDING TO YOUR WORKING ENVIRONMENT


########## Load Libraries ##########


library(tidyverse)
library(bbsBayes2)


########## Prepare the Script to Specific Species ##########


# Taking species name from the slurm script; if running locally, change spp_name

args = commandArgs(trailingOnly = TRUE)
spp_name = args[1]

min_max_route_years = 3 # default 3

voronoi_t = F # whether using Voronoi method to connect neighbors; default F

island_link_dist_factor = 1.2 # default 1.2; only applied when not using Voronoi method

# Print the arguments to verify

cat("Running model for species", spp_name, "with min_max_route_years =", min_max_route_years, "\n")

if (island_link_dist_factor != 1.2) {
  
  voronoi_t = F # not using Voronoi method
  cat("Finding neighbors using non-Voronoi method; 
      \nUsing island link distance factor =", island_link_dist_factor, "\n")
  
}

if (voronoi_t) {cat("Finding neighbors using Voronoi method \n")}

# Directories for saving outputs; change according to your working environment

output_dir = "Model_Fitting_Output/"
diagnos_dir = "Diagnostics/"
indices_dir = "Indices/"


########## Stratify Data ##########


s_grid = stratify(by = "latlong", species = spp_name)


########## Prepare Data ##########


# Spatial Models

p_grid = prepare_data(s_grid, min_n_routes = 1, min_max_route_years = min_max_route_years) # default 3
map = load_map(stratify_by = "latlong")
p_grid_spatial = prepare_spatial(p_grid, map, voronoi = voronoi_t, 
                                 island_link_dist_factor = island_link_dist_factor)

mod = prepare_model(p_grid_spatial, model = "first_diff", model_variant = "spatial",
                    calculate_log_lik = F)


########## Model Fitting ##########


# Default warmup and sampling iterations are 1000 each; for species that did not converge, double

fitted_model = run_model(mod, 
                         iter_warmup = 1000, 
                         iter_sampling = 1000,
                         output_dir = output_dir)

cat("Model for species", spp_name, "is completed\n")


########## Model Diagnostics ##########


summary_m = get_summary(fitted_model)

cat("Running diagnostics for species", spp_name, "...\n")

write_csv(summary_m, paste0(Diagnostics, gsub(" ", "", spp_name), ".csv"))

# Get the largest rhat

rhat_max = max(summary_m$rhat, na.rm = T)
rhat_max_n = max(summary_m %>% filter(str_detect(variable, "^n") & variable != "nu") %>% pull(rhat), 
                 na.rm = T)

if (rhat_max > 1.1) cat("Warning: maximum rhat for species", spp_name, "exceeds 1.1:", rhat_max, "\n")

cat("Diagnostics completed")


########## Extract Indices and Trends ##########


id_list = generate_indices(fitted_model)
indices = id_list[["indices"]]
indices$species = spp_name

cat("Writing relative abundance indices at\n",
          paste0(indices_dir, gsub(" ", "", spp_name), "_first_diff.csv\n"))

write_csv(indices, 
          paste0(indices_dir, gsub(" ", "", spp_name), "_first_diff.csv"))
