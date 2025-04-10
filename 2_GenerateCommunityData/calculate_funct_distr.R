# HELPER FUNCTION FOR CALCULATING CWM
# Hengxing Zou

calculate_funct_distr = function(tr, ab) {
  
  # if the abundance matrix has zero rows, return an empty data frame
  
  if (nrow(ab) == 0) {return(data.frame())}

  cwm = c()

  for (y in 1:nrow(ab)) {
    
    weights = ab[y, ] / sum(ab[y, ])
    
    # calculate cwm for year y
    
    cwm_year = weights %*% tr
    
    cwm = rbind(cwm, cwm_year)

  }
  
  cwm = as.data.frame(cwm)
  rownames(cwm) = rownames(ab)
  
  return(list(cwm))
  
}