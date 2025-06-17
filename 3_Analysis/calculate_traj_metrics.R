calculate_traj_metrics = function(trajs) {
  
  traj_char = foreach(r = names(trajs), 
                          .combine = rbind, .packages = c("tidyverse", "trajr")) %dopar% {
    
    mean_speed = mean(TrajDerivatives(trajs[[r]])$speed)
    sd_speed = sd(TrajDerivatives(trajs[[r]])$speed)
    mean_dc = mean(TrajDirectionalChange(trajs[[r]]))
    sd_dc = sd(TrajDirectionalChange(trajs[[r]]))
    Emax = TrajEmax(trajs[[r]])
    
    # Consider the start and end only to calculate effective displacement and angle
    
    x1 = trajs[[r]]$x[1]
    x2 = trajs[[r]]$x[51]
    y1 = trajs[[r]]$y[1]
    y2 = trajs[[r]]$y[51]
    
    eff_displ = sqrt((x2 - x1)^2 + (y2 - y1)^2)
    eff_angle = atan2(y2 - y1, x2 - x1)/pi*180
    angle_bins = case_when(eff_angle > 0 & eff_angle <= 90 ~ 1, 
                           eff_angle > 90 & eff_angle <= 180 ~ 2, 
                           eff_angle > -180 & eff_angle <= -90 ~ 3, 
                           eff_angle > -90 & eff_angle <= 0 ~ 4)
    
    # Output
    
    traj_char = tibble(region = r, 
                       mean_speed = mean_speed, 
                       sd_speed = sd_speed, 
                       mean_dc = mean_dc, 
                       sd_dc = sd_dc, 
                       Emax = Emax, 
                       eff_displ = eff_displ, 
                       eff_angle = eff_angle, 
                       angle_bins = angle_bins)
    
    traj_char
    
  }

}
