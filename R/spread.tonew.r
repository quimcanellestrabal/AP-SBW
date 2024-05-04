spread.tonew = function(land, nc, side, radius, outbreak, preoutbreak, w.wind, w.host,
                        reduc.nnew.outbreak, reduc.nnew.preoutbreak){
  
  # cat("SBW adjacent spreading", "\n" )
  
  ## Function to simulate spreading to cells not yet defoliated, that is, cells with ny.def0>=5 and tssbw>=30
  ## The future function sbw.spread.tonew(filter(land, ny.def0>=5)) will return a vector with cell.ids
  ## MB:  Prob.spread in cell c = Proportion host species in neighborhood x 
  ##      level of defoliation in neighborhood during previous year x some climatic variable 
  ##      (and perhaps wind direction in future versions)
  ## The question is, how many cells do I have to select from the pool of potential cells?
  ## It should be a number of cells proportional to the number of cells defoliated in the previous year¿?
browser()
  
  ## Compute SBW spreading potential according to A. species in the source cell,
  ## B. distance between the source and the target cells, and C. position of the target cells
  ## with respect to the main wind direction
  potential.spreading = sbw.spread.from.source(land, nc, radius)
  
  ## Select only those cells that at least one neighbor is defoliated
  potential = potential[potential$neigh.curr.def>0 & potential$x>0,]
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
  else if(nrow(potential)==1){
    sbw.new.sprd = potential$cell.id
  } else{
    sbw.new.sprd = integer()
  }
    
  return(sbw.new.sprd)
}
