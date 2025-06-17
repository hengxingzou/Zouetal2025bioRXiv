calculate_spatial_heterogeneity = function(positions, sf_map, comm_type = "base") {
  
  years = na.omit(unique(positions$year))
  
  positions_sf = right_join(positions, sf_map, by = c("region" = "ST_12"))
  
  sp_hetero = foreach(y = years, .packages = "tidyverse", .combine = rbind) %dopar% {
    
    sp_obj = positions_sf %>% 
      filter(year == y) %>% 
      filter(!is.na(quad)) %>% 
      sf::st_as_sf() %>% 
      sf::as_Spatial()
    
    coords = sp::coordinates(sp_obj)
    nb = spdep::tri2nb(coords, row.names = sp_obj@data[["region"]])
    nbw = spdep::nb2listw(nb, style = "W")
    gmoran = spdep::moran.test(sp_obj@data[["quad"]], nbw, alternative = "greater")
    
    as_tibble(as.list(gmoran[["estimate"]])) %>% 
      mutate(year = y, community = comm_type)
    
  }
  
  return(sp_hetero)
  
}
