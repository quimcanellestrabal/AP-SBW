spread.tonew = function(land, nc, side, radius, outbreak, preoutbreak, w.wind, w.host){
  
  # cat("SBW adjacent spreading", "\n" )
  
  ## Function to simulate spreading to cells not yet defoliated, that is, cells with ny.def0>=5 and tssbw>=30
  ## The future function sbw.spread.tonew(filter(land, ny.def0>=5)) will return a vector with cell.ids
  ## MB:  Prob.spread in cell c = Proportion host species in neighborhood x 
  ##      level of defoliation in neighborhood during previous year x some climatic variable 
  ##      (and perhaps wind direction in future versions)
  ## The question is, how many cells do I have to select from the pool of potential cells?
  ## It should be a number of cells proportional to the number of cells defoliated in the previous year¿?
  
  ## Compute neighborhood current defoliation and neighborhood host preference
  potential = neigh.influence.sbw.spread(land, nc, side, radius)
  
  ## The probability of sbw spread into a cell is proportional to neigh.curr.def, neigh.host.pref and 
  ## dominant wind directions
  ## Rescale the variables to the range [0,1] before applying any weight
  potential$x = w.wind * 0 + w.host * rescale(potential$neigh.host.pref, to=c(0,1)) +
    (1-pmax(w.wind+w.host,1)) * rescale(potential$neigh.curr.def, to=c(0,1))
  
  ## Select only those cells that at least one neighbor is defoliated
  potential = filter(potential, neigh.curr.def>0, x>0)
  if(nrow(potential)>1){
    if(outbreak<=10){  ## consolidated part of the outbreak
      #nnew = pmin(nrow(potential), sum(na.omit(land$curr.intens.def)>0)*runif(1, 0.15, 0.45)) #valors originals
      # nnew = pmin(nrow(potential), sum(na.omit(land$curr.intens.def)>0)*runif(1, 0.6, 1)) #valors meus perquè nnew sigui 10000-30000 de SAB i EPN
      num.new = nrow(potential) * (1-runif(1, 0, reduc.nnew.outbreak))
    }
    if(preoutbreak>0 | outbreak>10){  
      #nnew = pmin(nrow(potential), sum(na.omit(land$curr.intens.def)>0)*runif(1, 0.75, 1)) #valors originals
      # nnew = pmin(nrow(potential), sum(na.omit(land$curr.intens.def)>0)*runif(1, 1, 1.5))
      num.new = nrow(potential) * (1-runif(1, 0, reduc.nnew.preoutbreak))
    }
    sbw.new.sprd = sample(potential$cell.id, round(num.new), replace=F, prob=potential$x)
  }
  else
    sbw.new.sprd = potential$cell.id
  
  return(sbw.new.sprd)
}
