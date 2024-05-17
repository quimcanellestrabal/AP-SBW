########################################################################################################################
## Criteria that influence the SBW spreading to the neighbors:
## A. The species of the source cell
## B. The distance between source and target cells
## C. The position of the target cell with respect with the main wind direction
##########################################################################################

sbw.spread.from.source = function(land, nc, radius=12, side = 1){

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
  ## We use the weight_dist
  r = dist[,radius]
  weight_dist = round(exp(-dist^2/(r/2)^2),3)
  linear_dist = round((r-dist)/r,3)
  
  ## ---> Criteria A
  ## Matrix with the cell.id of the neighbor cells that are not defoliated yet
  neigh_nodefol = ids * (defol==0) 
  
  ## Copy in this data frame the cell.id of the source cells (1st column)
  neigh_nodefol[,1] = ids[,1]
  neigh_nodefol = data.frame(neigh_nodefol)
  names(neigh_nodefol)[1] = "source"  
  
  ## Relation between the source cell, the target cell andthe species of the source cell
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
  ## Élise wheel !  19/04/2024
  ## Do we need to turn the wheel according to a user-defined prevailing wind?
  ## We won't have a raster with the prevailing wind in each cell (we assume, too complex
  x <- c(radius:-radius)
  G <- expand.grid(x,x)
  names(G) <- c("x","y")
  # Cell id
  G$z <- nc*G$y + G$x #c(1:length(x)^2)
  # Distance from focal cell
  G$w <- sqrt((side*G$x)^2 + (side*G$y)^2)
  G$theta <- atan2(G$y, G$x) + 2*pi*abs(sign(G$y))*(1/2)*(1-sign(G$y))
  G$angle <- (180/pi)*G$theta
  # shift is the wind direction. 
  shift = 90 #From East to West
  G$shift <- (G$angle + shift) %% 360
  # if G$shift90 < 180 then G$w_angle = (180 - G$shift)/180
  # else G$w_angle = (360 - G$shift)/180
  G$w_angle <- abs(sign(G$shift-180))*((1/2)*(1-sign(G$shift-180))*((180 - G$shift)/180) + (1/2)*(1+sign(G$shift-180))*((360 - G$shift)/180))
  
  ## The position of each cell are in: G$z. And G$z is like position = ids[i,] - ids[i,1] (for the source cell 'i')
  source_wind = source_wspp %>% select(source, target) %>% mutate(position=target-source) %>% 
    left_join(select(G, z, w_angle), by=c("position"="z"))
    
  ## Merge the three criteria and compute the final weight 
  res = cbind(source_wspp, source_dist[,-1], source_wind[,"w_angle"]) %>% mutate(w=w_spp*w_dist*w_angle) %>%  
    group_by(target) %>% summarise(final_w=sum(w)) 
  
  ## Rescaling the final weight to [0,1]
  res$spread.potential = (res$final_w-min(res$final_w))/(max(res$final_w)-min(res$final_w)) 
  
  ## Questions:  
  ## Rescale the variables to the range [0,1] before applying any weight ??
  ## Do we want to apply a weight to each of these factors?
  ## Do we use a multiplicative or an additive formula?
  ## weight_factor_spp*w_spp + weight_factor_distance*w_dist + weight_factor_angel*w_angle
  
  return(res)  
}


