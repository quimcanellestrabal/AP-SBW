
#@ exemple:: load("data/landscape.rda"); land=landscape; preoutbreak=1; outbreak=1; calm=1; collapse=1; duration.last.outbreak=1

sbw.outbreak = function(land, params, tbls, preoutbreak=1, outbreak=1, calm=1, collapse=1, 
                        duration.last.outbreak=1, mask=NA){

  # 0.  Fix function  
  `%notin%` = Negate(`%in%`)

  ## Determine cell resolution in km
  cell.size.km = (mask@extent[2] - mask@extent[1])/ncol(mask) / 10^3
  
  ## Check that weights make sense
  if(params$w.wind>1){
    stop("Weight of wind factor cannot be greater than 1")
  }
  if(params$w.host>1){
    stop("Weight of neighbor host factor cannot be greater than 1")
  }
  
  ## 1. Defliation in the different phases of the SBW outbreak
  cat("   Defoliation in the ")
  if(preoutbreak>0){
    cat("pre-epidemic phase ", "\n")

    ## Potential are the cells són les cells en què hi pot haver sbw. Tota la funció s'aplica a "potential" i no a "land"
    potential = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"), temp>0.5, temp<2.8)
    
    # First epicenters
    epicenter = sample(potential$cell.id, size=rdunif(1,4,preoutbreak), replace=F,
                          prob=potential$ny.def0*(1200-potential$elev)/100)  # elevation threshold from Bouchard et al. 20xx
    sbw.new.sprd = epicenter
    
    ## Find between 20 to 40 neighs for teach epicenter and add to sbw.new.sprd.
    ## Each epicenter core has a different size indicated by the kmin and kmax parameters!
    for(i in 1:length(epicenter)){
      neighs = nn2(land[,c("x", "y")], land[land$cell.id==epicenter[i], c("x", "y")],
                   k=rdunif(1, params$kmin.bubble, params$kmax.bubble) , searchtype='priority')  
      nn.indx = neighs[[1]]
      sbw.new.sprd = unique(c(sbw.new.sprd, land$cell.id[nn.indx]))
    }

      # # just give some intensity to the core of the epicenters
      # land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
      #   sample(0:3, size=length(sbw.new.sprd), replace=T, prob=c(0.2,0.4,0.3,0.1))
      # 
      # # add some neighs of these epicenters #Quim radi a 8 (mida original = 12)
      # sbw.new.sprd = c(sbw.new.sprd, 
      #                  spread.tonew(land, nc=ncol(mask), side=cell.size.km, 
      #                               radius=8, outbreak, preoutbreak))
      # sbw.new.sprd = unique(sbw.new.sprd)
    
    # Select sbw.new.sprd only on potential cells
    potential = filter(land, spp %in% c("SAB", "EPN"))
    sbw.new.sprd = sbw.new.sprd[sbw.new.sprd %in% potential$cell.id]
    
    # and finally assign intensity to all of them (rewrite intensity just assigned to epicenter cores)
    land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
      sample(0:3, size=length(sbw.new.sprd), replace=T, prob=c(0.2,0.4,0.3,0.1))
  }

  if(outbreak>0){
    cat("epidemic phase ", "\n")
    ## Spatial spreading of the current outbreak to cells not yet defoliated, that is, cells with ny.def0>=5 & tssbw>=30
    ## The function 'spread.tonew' returns cell.ids
    ## The radius has to be variable to allow spreading further away or limit outbreak
    radius = rdunif(1, params$radius.outbreak.mid-params$radius.outbreak.range, params$radius.outbreak.mid+params$radius.outbreak.range) 
    sbw.new.sprd = spread.tonew(land, nc=ncol(mask), side=cell.size.km, radius=radius, outbreak, preoutbreak,
                                params$w.wind, params$w.host, params$reduc.nnew.outbreak, params$reduc.nnew.preoutbreak)
    sbw.new.sprd = unique(sbw.new.sprd)
    
    ## Only if some new cells are defoliated, assign level of defoliation
    if(length(sbw.new.sprd)>0){
      ## Select sbw.new.sprd only on potential cells
      # potential = filter(land, spp %in% c("SAB", "EPN"))
      sbw.new.sprd = sbw.new.sprd[sbw.new.sprd %in% potential$cell.id]
      
      ## Level of defoliation of the cells recently integrated in the outbreak (the sbw.new.spread cells)
      ## It can be 0 (no-defoliation), 1, 2 or 3!
      land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
        sample(0:3, size=length(sbw.new.sprd), replace=T, prob=c(0.2,0.4,0.3,0.1)) 
    }
  }

  if(collapse>0){ 
    cat("collapse phase ", "\n")
    ## if collapse, reduce number of new cells, only spontaneously
    potential = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"),temp>0.5, temp<2.8)
    sbw.new.sprd = sample(potential$cell.id, size=round(rlnorm(1, 2,1.5)), replace=F,  #mn(lnorm) ~ ifelse(collapse==1, 2, 1.5)
                          prob=potential$ny.def0*(1200-potential$elev)/100)
    
    # # #Some spontaneous new plots
    # if(rdunif(1,1,3)==1){
    #   radius = rdunif(1,6,12+collapse) #smaller radius than outbreak phase
    #   sbw.new.sprd = spread.tonew(land, nc=ncol(mask), side=2000/10^3,radius=radius, outbreak, preoutbreak)
    #   sbw.new.sprd = unique(sbw.new.sprd)
    #   #Select sbw.new.sprd only onn potential cells
    #   sbw.new.sprd = sbw.new.sprd[sbw.new.sprd %in% potential$cell.id]
    # }
    
    #Add intensity
    land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
      sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.5,0.3,0.2))
  }

  if(calm>0){    #calm>0 | collapse>0
    cat("calm phase ", "\n")
    ## not add new locations to the current outbreak if it is collapsing
    ## when calm, few spontaneously cells will have to be here and there defoliated
    potential = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"),temp>0.5, temp<2.8)
    sbw.new.sprd = sample(potential$cell.id, size=round(rlnorm(1, 2,1.5)), replace=F,  #mn(lnorm) ~ ifelse(collapse==1, 2, 1.5)
                          prob=potential$ny.def0*(1200-potential$elev)/100)
    #Add intensity
    land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
      sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.5,0.3,0.2))
  }
  
  
  ## 2. Update SBW tracking variables for newly defoliated cells
  land$ny.def[land$cell.id %in% sbw.new.sprd] = 
    ifelse(land$curr.intens.def[land$cell.id %in% sbw.new.sprd]==0, 0, 1)
  land$ny.def0[land$cell.id %in% sbw.new.sprd] = 
    ifelse(land$curr.intens.def[land$cell.id %in% sbw.new.sprd]>0, 0, 
           land$ny.def0[land$cell.id %in% sbw.new.sprd]+1)
  land$cum.intens.def[land$cell.id %in% sbw.new.sprd] = 
    land$curr.intens.def[land$cell.id %in% sbw.new.sprd]
  
  
  ## 3. Update level of defoliation intensity of already defoliated cells, that is, 
  ## cells with ny.def0<5 & tssbw>=30, and ny.def<=18 --> otherwise, the defoliation is forced to stop.
  sbw.curr.sprd = land$cell.id[land$ny.def0<5 & land$tssbw>=30 & land$ny.def<=18 & 
                                 land$cell.id %notin% sbw.new.sprd]
  curr.outbreak = filter(land, cell.id %in% sbw.curr.sprd)
  ## Level of defoliation of these cells. It can be 0 (no-defol), 1, 2 or 3!
  land$curr.intens.def[land$cell.id %in% sbw.curr.sprd] = 
    intens.def.curr(filter(land, cell.id %in% sbw.curr.sprd), params, tbls, preoutbreak, outbreak, collapse, calm)
  ## Update SBW tracking variables
  land$ny.def[land$cell.id %in% sbw.curr.sprd] = land$ny.def[land$cell.id %in% sbw.curr.sprd]+
    ifelse(land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]==0, 0, 1)
  land$ny.def0[land$cell.id %in% sbw.curr.sprd] = 
    ifelse(land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]>0, 0, 
           land$ny.def0[land$cell.id %in% sbw.curr.sprd]+1)
  land$cum.intens.def[land$cell.id %in% sbw.curr.sprd] = 
    land$cum.intens.def[land$cell.id %in% sbw.curr.sprd]+land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]  
  
  
  ## 4. Stop defoilating cells that have been defoliated for more than 18 years
  ## Warning, not consider cells that have just been defoliated
  ## And avoid recurrence of defoliation in this oubreak by making tssbw==0
  sbw.stop = land$cell.id[land$ny.def0==0 & land$ny.def>18 & land$cell.id %notin% sbw.curr.sprd]
  land$curr.intens.def[land$cell.id %in% sbw.stop] = 0
  land$ny.def0[land$cell.id %in% sbw.stop] = 1
  land$tssbw[land$cell.id %in% sbw.stop] = 0
  
  
  ## 5. Increase number of year of non-defoliation for the remaing cells
  land$ny.def0[land$cell.id %notin% c(sbw.new.sprd, sbw.curr.sprd, sbw.stop)] = 
    land$ny.def0[land$cell.id %notin% c(sbw.new.sprd, sbw.curr.sprd, sbw.stop)]+1
  
  
  ## 6. For those cells that have just accumulated 5 years of non-defoliation, 
  ## reset the counter of number of years of defoliation and 
  ## cumulative intensity of defoliaion
  land$ny.def[land$ny.def0==5] = 0
  land$cum.intens.def[land$ny.def0==5] = 0
  
  
  ## 7. Kill cells according to number of years of defoliation and species composition
  ## with probability as cumulative*current intensity of defoliation
  if(outbreak>0 | collapse>0){
    kill.cells = forest.mortality(land)
  }
  else{
    kill.cells = integer()
  }
  
  
  ## 8. Set outbreak parameters
  if(preoutbreak>0){ 
    preoutbreak = preoutbreak-1
    phase = "preoutbreak"
    if(preoutbreak==0)
      duration.last.outbreak = outbreak = rdunif(1,0,2)+params$duration.last.outbreak  # +6 Gray2008 #Quim: ho he canviat (outbreak entre 9 i 11), l'original era "rdunif(1,-1,2) (outbreak entre 8 i 10"
  }
  else if(outbreak>0){ 
    outbreak = outbreak-1
    phase = "outbreak"
    if(outbreak==0)
      collapse = rdunif(1,4,6) #Original entre 3 i 6
  }
  else if(collapse>0){  
    phase = "collapse"
    collapse = collapse-1
    if(collapse==0) 
      calm = round(rnorm(1, 15, 2))
  }
  else if(calm>0){
    phase = "calm"
    calm = calm-1    
    if(calm==0)
      preoutbreak = rdunif(1,3,4)
  }
  
  res = list(kill.cells=kill.cells, land.sbw=land, preoutbreak=preoutbreak, outbreak=outbreak, collapse=collapse, calm=calm)
  return(res) 
} 

