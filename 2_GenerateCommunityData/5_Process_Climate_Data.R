library(tidyverse)
library(foreach)
library(doSNOW)

cl = makeCluster(10, outfile = "")
registerDoSNOW(cl)

# Change the directory according to your working environment

database_dir = "Database/"


########## Calculate Bioclimatic Variables ##########


# Read data from WorldClim2
# Change this to appropriate file path where the data is hosted

precip_df = read_csv("MonthlyPrecip.csv")
tmax_df = read_csv("MonthlyTmax.csv")
tmin_df = read_csv("MonthlyTmin.csv")

all_bioclim = foreach(g = unique(precip_df$ST_12), 
                      .combine = rbind, .packages = "tidyverse") %dopar% {
  
  output_g = list()
  
  for (y in unique(precip_df$year)) {
          
    cat("Generating bioclimatic variables for", g, "in", y, "\n")
    
    precip_y_g = precip_df %>% 
      filter(year == y, ST_12 == g)
    
    tmax_y_g = tmax_df %>% 
      filter(year == y, ST_12 == g)
    
    tmin_y_g = tmin_df %>% 
      filter(year == y, ST_12 == g)
    
    bioclim = as_tibble(dismo::biovars(precip_y_g$mean_data, tmin_y_g$mean_data, tmax_y_g$mean_data)) %>% 
      mutate(year = y, ST_12 = g)
      
    output_g[[y]] = bioclim
    
  }
  
  bind_rows(output_g)
  
}

# Write data

write_csv(all_bioclim, paste0(database_dir, "Bioclim.csv"))
  
