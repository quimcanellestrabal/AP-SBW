neigh.influence.sbw.spread = function(land, nc=804, side=2, radius=12){
  
  # radius = 1
  # y = c(-radius:radius); y
  # z = c(y, ncol(MASK)+y, -(ncol(MASK)+y)); z[z!=0]
  # 
  # radius = 2
  # y = c(-radius:radius); y
  # z = c(y, ncol(MASK)+y, -(ncol(MASK)+y), 2*ncol(MASK)+y, -(2*ncol(MASK)+y)); z[z!=0]
  
  ## Neighbor position and distance to the focal cell in km
  # x = c(-radius:radius)
  # z = x; z
  # w = abs(z)*side; w
  # for(i in 1:radius){
  #   z = c(z, i*nc+x, -(i*nc+x)) 
  #   w = c(w, sqrt((abs(x)*side)^2+(i*side)^2), sqrt((abs(x)*side)^2+(i*side)^2) )
  # }
  # z = z[z!=0]
  # w = w[w!=0]
  
  x <- c(radius:-radius)
  G <- expand.grid(x,x)
  names(G) <- c("x","y")
  
  # Cell id
  G$z <- nc*G$y + G$x #c(1:length(x)^2)
  
  # Distance from focal cell
  G$w <- sqrt((side*G$x)^2 + (side*G$y)^2)
  
  # Angles
  # angles in the neighborhood are assumed to correspond to E=0, N=90, W=180, and S=270
  # atan2(y,x) is a type of atan that returns 0 even if x is 0 (no NaN value)
  # for cells in the 1st and 2nd quadrants it returns angles between 0 and pi, 
  # whereas for cells in the 3rd and 4th quadrants it returns angles between 0 and -pi
  # To have angles between 0 and 360 you need to
  # 1) convert radians to degrees: multiply by (180/pi)
  # 2) Add 360 (2*pi) to angles for cells in the 3rd and 4th: 
  # i.e. only for cells with negative value of y: the term (1/2)*(1-sign(G$y)) ensures that (when y>0, this term is 0, otherwise it is 1)
  # and only for cells with strictly negative value of y (i.e. not y=0): the term abs(sign(G$y)) ensures that (when y=0, this term is 0)
  
  G$teta <- atan2(G$y, G$x) + 2*pi*abs(sign(G$y))*(1/2)*(1-sign(G$y))
  G$teta <- (180/pi)*G$teta
  
  # Remove focal cell (z = 0)
  G <- filter(G, z!=0)
  
  ## Look at current level of the defoliation in the neighborhood of the cells not 
  ## currently defoliated and that last mortality by outbreak is at least 30 years.
  potential = filter(land, ny.def0>=5, tssbw>=30)
  nslice = 10
  upper = round(nrow(potential)/nslice)
  ## First slice
  potential.slice = potential[1:upper,]
  neigh.curr.def = .compute.neigh.curr.def(land, potential.slice, G$z, G$w)
  neigh.host.pref = .compute.neigh.host.pref(land, potential.slice, G$z,G$w)
  ## Second to n-1 slice
  for(i in 1:(nslice-1)){
    # cat(i, "\n")
    potential.slice = potential[(i*upper+1):((i+1)*upper),]
    neigh.curr.def = rbind(neigh.curr.def, .compute.neigh.curr.def(land, potential.slice, G$z, G$w))
    neigh.host.pref = rbind(neigh.host.pref, .compute.neigh.host.pref(land, potential.slice, G$z, G$w))
  }
  ## Last slice
  # cat("last")
  potential.slice = potential[((nslice-1)*upper+1):nrow(potential),]
  neigh.curr.def = rbind(neigh.curr.def, .compute.neigh.curr.def(land, potential.slice, G$z, G$w))
  neigh.host.pref = rbind(neigh.host.pref, .compute.neigh.host.pref(land, potential.slice, G$z, G$w))
  
  ## Aggregate all the info
  dta = data.frame(neigh.curr.def, neigh.host.pref$x)
  names(dta) = c("cell.id", "neigh.curr.def", "neigh.host.pref")
  
  return(dta)
}


.compute.neigh.curr.def = function(land, potential.slice, z, w){
  nneigh = length(z)
  neighs = data.frame(cell.id=rep(potential.slice$cell.id, each=nneigh),
                       neigh.id=rep(potential.slice$cell.id, each=nneigh) + rep(z, nrow(potential.slice)),
                       w=rep(w, nrow(potential.slice))) %>% 
    filter(neigh.id %in% land$cell.id)
  neigh.def = left_join(neighs, dplyr::select(land, cell.id, curr.intens.def), by=c("neigh.id"="cell.id")) %>% 
    group_by(cell.id) %>% summarise(x=sum(curr.intens.def/(w*0.5)))
  return(neigh.def)
}

.compute.neigh.host.pref = function(land, potential.slice, z, w){
  nneigh = length(z)
  breaks = c(0,20,40,60,80,100,999)
  tags = c("C10", "C30", "C50", "C70", "C90", "OLD")
  neighs = data.frame(cell.id=rep(potential.slice$cell.id, each=nneigh),
                       neigh.id=rep(potential.slice$cell.id, each=nneigh) + rep(z, nrow(potential.slice)),
                       w=rep(w, nrow(potential.slice))) %>% 
    filter(neigh.id %in% land$cell.id)
  neigh.host = left_join(neighs, dplyr::select(land, cell.id, spp, age), by=c("neigh.id"="cell.id")) %>% 
    mutate(age.class=cut(age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)) %>% 
    mutate(host.pref=ifelse(spp=="SAB" & age.class %in% c("C50", "C70", "C90", "OLD"), 1,
                            ifelse(spp=="SAB" & age.class %in% c("C10", "C30"), 0.75,
                                   ifelse(spp=="EPN" & age.class %in% c("C50", "C70", "C90", "OLD"), 0.5,
                                          ifelse(spp=="EPN" & age.class %in% c("C10", "C30"), 0.25, 0)))) ) %>% 
    group_by(cell.id) %>% summarise(x=sum(host.pref/(w*0.5)))
  return(neigh.host)
}
