library(tidyverse)
library(raster)
library(terra)
library(exactextractr)
library(sf)

library(foreach)
library(doParallel)
# Not using doSNOW like other scripts here because it creates a consistent error

# cl = makeCluster(10, outfile = "")
registerDoParallel(cores = 10)

# Change the directory according to your working environment

database_dir = "Database/"

# Read all anthromes (.asc raster files)
# Change this to appropriate file path where the data is hosted

files = Sys.glob("../Databases/HYDE/Anthromes/*.asc")

all_rasters = lapply(files, function(x) raster(x))

names(all_rasters) = 1970:2021

# Read BBS spatial extent

BBS = st_read(paste0(database_dir, "BBS_LatLong_strata.shp"))

BBS_reproj = st_transform(BBS, crs = crs(all_rasters[[1]])) %>% # this is the crs of all anthromes data
  mutate(Grid_ID = 1:nrow(.))

BBS_ext = st_bbox(BBS_reproj)

# Crop the original raster for faster processing
  
cropped_all_rasters = lapply(all_rasters, function(x) terra::crop(x, BBS_ext))

# All categories of anthromes in North America; remove NAs
# No change of categories over the years
# The only category missing from the global anthromes is 21 (Village - Rice)

all_cat = unique(cropped_all_rasters$`1970`@data@values)[-1]

# For each year, aggregate anthromes to grid level, then calculate coverage of each anthrome in a grid

all_extractions = foreach(y = names(all_rasters), .packages = "tidyverse", .combine = rbind) %dopar% {
  
  cropped_r = cropped_all_rasters[[y]]
  
  extraction = exactextractr::exact_extract(cropped_r, BBS_reproj, 
                                            fun = function(values, coverage_fraction) {
    
    # Calculate the proportion of each category
    coverage = numeric() # empty vector to hold calculated proportions
    
    for (c in all_cat) {
      
      # Remove all NA in sums
      coverage[as.character(c)] = sum(coverage_fraction[values == c], na.rm = T) / sum(coverage_fraction, na.rm = T)
      
    }
    
    return(coverage)
    
  }, coverage_area = T, progress = T) %>% t() %>% data.frame()
  
  extraction %>% 
    mutate(Year = as.numeric(y)) %>% 
    rownames_to_column("Grid_ID") %>% 
    mutate(Grid_ID = as.numeric(Grid_ID))
    
}

# Match grids to actual coordinates

all_extractions_sf = all_extractions %>% 
  left_join(BBS_reproj, by = "Grid_ID") %>% 
  pivot_longer(starts_with("X"), names_to = "Anthrome", values_to = "Coverage") %>%
  mutate(Anthrome = str_sub(Anthrome, 2, nchar(Anthrome)))

# Write extracted proportions for all grids and years

write_csv(st_drop_geometry(all_extractions_sf) %>% 
            dplyr::select(-geometry), 
          paste0(database_dir, "Anthromes_Grid_Years_Prop.csv"))

# Compare with original raster to cross check

# plot(cropped_all_rasters$`1970`)

all_extractions_sf %>%
  filter(Year == 1970) %>% 
  filter(Anthrome %in% c("31", "32", "33", "34")) %>%
  group_by(Grid_ID) %>%
  mutate(Coverage = sum(Coverage)) %>%
  ggplot(aes(fill = Coverage, geometry = geometry)) +
  geom_sf()
