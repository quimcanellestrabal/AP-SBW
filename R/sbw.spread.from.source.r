########################################################################################################################
## Criteria that influence the SBW spreading to the neighbors:
## A. The species of the source cell
## B. The distance between source and target cells
## C. The position of the target cell with respect with the main wind direction
##########################################################################################


sbw.spread.from.source = function(land, nc, wind_dir = 90, radius=12, side = 1){  ##From Est (works for different angle between 0 and 360


  ## Select all the cells that are currently defoliated  
  source.cells = land[land$ny.def>0,]
  
  ## Find the neighbours of all source cells in a radius
  neighs = nn2(land[,c("x", "y")], land[land$cell.id %in% source.cells$cell.id, c("x", "y")],
               k=radius, searchtype='priority')  
  nn.indx = neighs[[1]] # position of cells within the dataframe 'land'
  nn.dist = neighs[[2]] # distance between source and target of cells
  
  ## We need to find out 
  ## 1. the species of the source cell associated to each neighbor (spps)
  ## 2. if the neighbors are defoliated or not (defol)
  ## 3. the distance between the target and the source cells (dist)
  ids = matrix(land$cell.id[nn.indx], nrow=nrow(source.cells), ncol=radius)  # the real id of each cell (it is different from its position within the 'land' dataframe)
  spps = matrix(land$spp[nn.indx], nrow=nrow(source.cells), ncol=radius)
  defol = matrix(land$ny.def[nn.indx], nrow=nrow(source.cells), ncol=radius)
  dist = matrix(nn.dist, nrow=nrow(source.cells), ncol=radius)
  
  ## Maximum radius for each target cell and the corresponding weighted and linear distances
  ## We use the weight_dist  # EF: NOT USED. REMOVE?
  r = dist[,radius]
  weight_dist = round(exp(-dist^2/(r/2)^2),3)
  linear_dist = round((r-dist)/r,3)
  
  ## ---> Criteria A
  ## Matrix with the cell.id of the neighbor cells that are not defoliated yet
  neigh_nodefol = ids * (defol==0) 
  
  ## Copy in this data frame the cell.id of the source cells (1st column)
  neigh_nodefol[,1] = ids[,1]
  neigh_nodefol = data.frame(neigh_nodefol)
  names(neigh_nodefol)[1] = "source"  ## EF: AND THE OTHER ARE THE TARGET CELLS?
  
  ## Relation between the source cell, the target cell and the species of the source cell
  source_wspp = pivot_longer(neigh_nodefol, cols=2:(radius), values_to="target") %>% select(-name) %>% 
    filter(target!=0) %>% left_join(select(land, cell.id, spp), by=c("source" = "cell.id")) %>% 
    rename(source_spp = spp)
  
  ## Apply the weight corresponding to these species
  source_wspp$w_spp = ifelse(source_wspp$source_spp %in% c("EPN", "SAB"), 0.7, 0.3)
  
  ## ---> Criteria B
  ## Matrix with the distance of the neighbor cells that are not defoliated yet
  neigh_dist = dist * (defol==0) 
  
  ## Copy in this data frame the cell.id of the source cells (1st column)
  neigh_dist[,1] =  ids[,1]
  neigh_dist = data.frame(neigh_dist)
  names(neigh_dist)[1] = "source" 
  
  ## Relation between the source cell and the distance with the target cells (not explicitly indicated in this data frame,
  ## but we keep the order of them)
  source_dist =  pivot_longer(neigh_dist, cols=2:(radius), values_to="dist") %>% select(-name)  %>% filter(dist!=0)
  r = max(source_dist$dist)-2000
  
  ## Apply the weight corresponding to the distance
  source_dist$w_dist = exp(-(source_dist$dist-2000)^2/(r/2)^2)
  
  ## ---> Criteria C  
  ## Élise updated 14-05-2024
  
  x <- c(radius:-radius)
  G <- expand.grid(x,x)
  names(G) <- c("x","y")
  # Cell id
  G$z <- nc*G$y + G$x #c(1:length(x)^2)
  # Distance from focal cell
  G$w <- sqrt((side*G$x)^2 + (side*G$y)^2)
  # Each cell of the neighborhood is given an angle according to the convention: N=0, E=90, S=180, W=270
  G$theta <- atan2(G$x, G$y) + 2*pi*abs(sign(G$x))*(1/2)*(1-sign(G$x)) #in rads
  G$angle <- (180/pi)*G$theta #in degrees
  
  # Weight wind: 
  # cells in the direction of the prevailing wind (wind_dir) are given higher weights
  # Presenting wind_dir is a user-define number identical for all cells.
  # However, it could be a raster with the prevailing wind in each cell (verification that the code still works would be needed)
  # Weights decrease on cells away from the preferred wind direction

  G$wind_left <- (G$angle - wind_dir) %% 360
  G$wind_right <- (wind_dir - G$angle) %% 360
  G$wind <- G$wind_left*(G$wind_left<=180) + G$wind_right*(G$wind_right<180)
  G$w_wind <- (1/180)*G$wind  # test: matrix(G$w_wind, nrow = 5, byrow = TRUE)
  # Remove focal cell as w_wind is not correct at the center
  G <- filter(G, z!=0)
  

  ## The position of each cell are in: G$z. And G$z is like position = ids[i,] - ids[i,1] (for the source cell 'i')
  source_wind = source_wspp %>% select(source, target) %>% mutate(position=target-source) %>% 
    left_join(select(G, z, w_wind), by=c("position"="z"))

  ## Merge the three criteria and compute the final weight 
  res = cbind(source_wspp, source_dist[,-1], source_wind[,"w_wind"]) %>% 
    mutate(w_add = 1/3*w_spp + 1/3*w_dist + 1/3*w_wind, w_multi=w_spp*w_dist*w_wind) %>%  
    group_by(target) %>% summarise(final_w_add=sum(w_add), final_w_multi=sum(w_multi)) 

  ## Rescaling the final weight to [0,1]
  res$spread_potential_add = (res$final_w_add-min(res$final_w_add))/(max(res$final_w_add)-min(res$final_w_add)) 
  res$spread_potential_multi = (res$final_w_multi-min(res$final_w_multi))/(max(res$final_w_multi)-min(res$final_w_multi)) 
  
  ## Questions:  
  ## Rescale the variables to the range [0,1] before applying any weight ??
  ## Do we want to apply a weight to each of these factors?
  ## Do we use a multiplicative or an additive formula?
  ## weight_factor_spp*w_spp + weight_factor_distance*w_dist + weight_factor_angel*w_wind
  
  return(res)  
}


