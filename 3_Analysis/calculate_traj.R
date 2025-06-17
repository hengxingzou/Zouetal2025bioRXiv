calculate_traj = function(pca_grids) {
  
  grids = unique(pca_grids$region)
  
  trajs = foreach(r = grids, .packages = c("tidyverse", "trajr")) %dopar% {
    
    coord = pca_grids %>% 
      filter(region == r)
    
    traj = TrajFromCoords(coord, xCol = 4, yCol = 5, fps = 1)
    
    traj
    
  }
  
  names(trajs) = grids
  return(trajs)
  
}

