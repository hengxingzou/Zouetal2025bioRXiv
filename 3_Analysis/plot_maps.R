plot_maps = function(endpoints, traj_char, comm_type = "base") {
  
  comm_label = case_when(comm_type == "base" ~ "Complete", 
                         comm_type == "topspp" ~ "Top", 
                         comm_type == "restspp" ~ "Standard")
  
  p_start_end_quad =
    ggplot() + 
    geom_sf(data = endpoints, 
            aes(fill = quad, geometry = geometry), 
            color = "gray", size = 0.1, 
            inherit.aes = F) + 
    coord_sf() + 
    scale_fill_gradientn(colors = color_4) +
    ggtitle(paste0(comm_label, ", Quadrants of Endpoints")) +
    facet_wrap(year ~ ., ncol = 1) + 
    theme_bw() + 
    theme(title = element_text(size = 12), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.position = "none")
  
  p_mean_speed_map =
    ggplot() + 
    geom_sf(data = traj_char, 
            aes(fill = sqrt(mean_speed), geometry = geometry), 
            color = "gray", size = 0.1, 
            inherit.aes = F) + 
    coord_sf() + 
    scale_fill_gradientn(colors = color_4, 
                         limits = range_mean_speed) +
    ggtitle(paste0(comm_label, ", Mean Speed")) +
    theme_bw() + 
    theme(title = element_text(size = 12), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.position = "none")
  
  p_eff_angle_map =
    ggplot() + 
    geom_sf(data = traj_char, 
            aes(fill = as.factor(angle_bins), geometry = geometry), 
            color = "gray", size = 0.1, 
            inherit.aes = F) + 
    coord_sf() + 
    scale_fill_manual(name = "Effective Angle", values = color_4) +
    ggtitle(paste0(comm_label, ", Effective Angle")) +
    theme_bw() + 
    theme(title = element_text(size = 12), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.position = "none")
  
  p_Emax_map =
    ggplot() + 
    geom_sf(data = traj_char, 
            aes(fill = Emax, geometry = geometry), 
            color = "gray", size = 0.1, 
            inherit.aes = F) + 
    coord_sf() + 
    scale_fill_gradientn(colors = color_4, 
                         limits = range_Emax) +
    ggtitle(paste0(comm_label, ", Directionality (Emax)")) +
    theme_bw() + 
    theme(title = element_text(size = 12), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.position = "none")
  
  p_eff_displ_map =
    ggplot() + 
    geom_sf(data = traj_char, 
            aes(fill = sqrt(eff_displ), geometry = geometry), 
            color = "gray", size = 0.1, 
            inherit.aes = F) + 
    coord_sf() + 
    scale_fill_gradientn(colors = color_4, 
                         limits = range_eff_displ) +
    ggtitle(paste0(comm_label, ", Effective Displacement")) +
    theme_bw() + 
    theme(title = element_text(size = 12), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.position = "none")

  p_start_end_angle =     
    ggplot() + 
    geom_sf(data = endpoints %>% 
              select(year, region, quad, LONGD:geometry) %>% 
              rbind(traj_char %>% 
                     mutate(year = "Effective Angle") %>% 
                     rename(quad = angle_bins) %>% 
                     select(year, region, quad, LONGD:geometry)),
            aes(fill = quad, geometry = geometry), 
            color = "gray", size = 0.1, 
            inherit.aes = F) + 
    coord_sf() + 
    scale_fill_gradientn(colors = color_4) +
    # ggtitle(paste0(comm_label)) +
    facet_wrap(year ~ ., ncol = 1, ) + 
    theme_bw() + 
    theme(title = element_text(size = 12), 
          panel.grid = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text = element_blank(), 
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 10), 
          legend.position = "none", 
          strip.background = element_blank(), 
          strip.text = element_blank())
  
  pls = list(p_start_end_quad, p_mean_speed_map, p_eff_angle_map, p_Emax_map, p_eff_displ_map, p_start_end_angle)
  
  names(pls) = c("start_end", "mean_speed", "eff_angle", "Emax", "eff_displ", "start_end_angle")
  
  return(pls)
  
}