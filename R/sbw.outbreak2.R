
#@ exemple:: load("data/landscape.rda")
#@ exemple:: land=landscape; preoutbreak=1; outbreak=1; calm=1; collapse=1; duration.last.outbreak=1

sbw.outbreak = function(land, preoutbreak=1, outbreak=1, calm=1, collapse=1, duration.last.outbreak=1){

  
      cat("SBW outbreak: ", "\n")
      
      ## 1. Defliation in the different phases of the SBW outbreak
      cat("Defoliation in the ")
      if(preoutbreak>0){
        cat("pre-epidemic phase ", "\n")
        potential = filter(land, ny.def0>=5, tssbw>=30, temp>0.5, temp<2.8)
        # preoutbreak=5
        sbw.new.sprd = sample(potential$cell.id, size=rdunif(1,1,6-preoutbreak), replace=F,
                               prob=potential$ny.def0*(1200-potential$elev)/100)
        ## find between 20 to 40 neighs and add to sbw.new.sprd
        neighs = nn2(dplyr::select(land,x,y), filter(land, cell.id %in% sbw.new.sprd) %>% dplyr::select(x,y),
                  k=rdunif(1,20,40) , searchtype='priority')
        # neighs = nn2(dplyr::select(land,x,y), filter(land, cell.id %in% sbw.new.sprd) %>% dplyr::select(x,y),
        #               k=40 , searchtype='priority')  #test
        nn.indx = neighs[[1]]
        sbw.new.sprd = land$cell.id[nn.indx]
        # just give some intensity to the core of the epicenters
        land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
          sample(0:3, size=length(sbw.new.sprd), replace=T, prob=c(0.2,0.4,0.3,0.1))
        # add some neighs of these epicenters
        sbw.new.sprd = c(sbw.new.sprd, 
            spread.tonew(land, nc=ncol(mask), side=res(mask)[1]/10^3, 
                             radius=12, outbreak, preoutbreak))
        sbw.new.sprd = unique(sbw.new.sprd)
        # and finally assign intensity to all of them (rewrite intensity just assigned to epicenter cores)
        land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
          sample(0:3, size=length(sbw.new.sprd), replace=T, prob=c(0.2,0.4,0.3,0.1))
      }
      # else 
      if(outbreak>0){
        cat("epidemic phase ", "\n")
        ## Spatial spreading of the current outbreak to cells not yet defoliated, that is,
        ## cells with ny.def0>=5 & tssbw>=30
        ## The function 'spread.tonew' returns cell.ids.
        ## once in a while vary the radius, to allow spreading furhter away, or reduce it to limit outbreak....
        radius = rdunif(1,2,15) # 4 to 60 km
        sbw.new.sprd = spread.tonew(land, nc=ncol(mask), side=2000/10^3,    #**en comptes de 2000 podriem posar raster::res(mask)[1], però el codi fa el tonto!
                                         radius=radius, outbreak, preoutbreak)
        sbw.new.sprd = unique(sbw.new.sprd)
        ## Level of defoliation of the cells recently integrated in the outbreak (the sbw.new.spread cells)
        ## It can be 0 (no-defol), 1, 2 or 3!
        land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
          sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.1,0.8,0.1))
      }
      #else 
      if(calm>0 | collapse==1){ 
        cat("calm phase ", "\n")
        ## not add new locations to the current outbreak if it is collapsing
        ## when calm, few spontaneously cells will have to be here and there defoliated
        potential = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"), temp>0.5, temp<2.8)
        sbw.new.sprd = sample(potential$cell.id, size=round(rlnorm(1, 2,1.5)), replace=F,  #mn(lnorm) ~ ifelse(collapse==1, 2, 1.5)
                               prob=potential$ny.def0*(1200-potential$elev)/100)
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
        
      ## 8. Set outbreak parameters
      if(preoutbreak>0 & done){ # pre-epidemic
        preoutbreak = preoutbreak-1
        phase = "preoutbreak"
        if(preoutbreak==0)
          duration.last.outbreak = outbreak = rdunif(1,-1,1)+params$duration.last.outbreak  # +6 Gray2008
        done = FALSE
      }
      else if(outbreak>0 & done){ # epidemia
        outbreak = outbreak-1
        phase = "outbreak"
        if(outbreak==0)
          collapse = rdunif(1,3,5)
        done = FALSE
      }
      else if(collapse>0 & done){  # collapse
        phase = "collapse"
        collapse = collapse-1
        if(collapse==0) #finishing the collapse
          calm = round(rnorm(1, 15, 1))
        done = FALSE
      }
      else if(calm>0 & done){
        phase = "calm"
        calm = calm-1    
        if(calm==0)
          preoutbreak = rdunif(1,3,4)
      }
      # cat(phase, "\n")
      done = TRUE
      

      res = list(kill.cells=kill.cells, land.sbw=land, preoutbreak=preoutbreak, outbreak=outbreak, collapse=collapse, calm=calm)
      return(res) 
      
      
  #     return(kill.cells)  
  # 
  # assign("land.sbw", land, envir=.GlobalEnv)
  # assign("preoutbreak", preoutbreak, envir=.GlobalEnv)
  # assign("outbreak", outbreak, envir=.GlobalEnv)
  # assign("collapse", collapse, envir=.GlobalEnv)
  # assign("calm", calm, envir=.GlobalEnv)
      
} 

