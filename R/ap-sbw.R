#' Quebec Landscape Dynamic Model: Artificial planting and SBW
#'
#' Run the landscape dynamic model QLDM that includes the processes of spruce budworm outbreaks, clear-cuts,
#' artificial planting, post-disturbance regeneration, and natural forest succession.
#'
#' @param scn Indicates the scn (setted in default.tables) that defines harvesting and artificial planting rates. 'scn0':'scn8'
#' @param is.sbw A flag to indicate that spruce budworm outbreaks are a simulated process of the model
#' @param is.harvesting A flag to indicate that harvesting by clear and partial cuts is a simulated process of the model
#' @param custom.params List with the model paramaters and default and/or user-defined values 
#' @param rcp Climate projection, either \code{NA} (default), 'rcp45' or 'rcp85' 
#' @param nrun Number of replicates to run the model
#' @param time.step Number of years of each time step
#' @param time.horizon Number of years of the model simulation, it has to be a multiple \code{time.step}
#' @param save.land A flag to save as a RDS file the \code{landscape} data frame at the time step indicated in \code{out.seq}
#' @param time.save Numeric vector with the time steps the \code{landscape} is saved
#' @param out.path String with the directory path to save the \code{landscape} data frame at each time step indicated in \code{out.seq}
#'
#' @return A list with the following items:
#'  \itemize{
#'    \item{\code{SppByAgeClass}: A data frame of species abundance per age class, with columns:
#'      \itemize{
#'         \item{\code{run}: Number of replicate.}
#'         \item{\code{year}: Year YYYY.}
#'         \item{\code{mgmt.unit}: Code of the forest management unit (FMU).}
#'         \item{\code{spp}: Code of the species.}
#'         \item{\code{age.class}: Code of the age class, \code{C10}, \code{C30}, \code{C50}, \code{C70}, \code{C90}, and \code{OLD}.}
#'         \item{\code{area}: Area in \eqn{km^{2}}.}
#'       }
#'    }
#'    \item{\code{SuitabilityClass}: A data frame of suitability of potential species per bioclimatic domain, with columns:
#'      \itemize{
#'         \item{\code{run}: Number of replicate.}
#'         \item{\code{year}: Year YYYY.}
#'         \item{\code{bioclim.domain}: Code of the bioclimatic domain.}
#'         \item{\code{spp}: Code of the species.}
#'         \item{\code{poor}: Area in \eqn{km^{2}} of poor environmental suitability.}
#'         \item{\code{med}: Area in \eqn{km^{2}} of intermediate environmental suitability.}
#'         \item{\code{good}: Area in \eqn{km^{2}} of good environmental suitability.}
#'       }
#'    }
#'  }
#'  
#'  
#'@example:: scn="scn1"; nrun=3; time.step=1; time.horizon=20; save.land=T; time.save=5; out.seq=NA; out.path=paste0("outputs/test2"); custom.params=NA; rcp='rcp85'; is.sbw=T; is.harvesting=T
#'@example:: scn="scn1"; is.sbw=T; is.harvesting=T; custom.params=NA; rcp='rcp85'; nrun=3; time.step=1; time.horizon=20; save.land=F; time.save=NA; out.path=NA 
#'@example:: kk=ap.sbw(scn="scn1", is.sbw=T, is.harvesting=T, custom.params=NA, rcp='rcp85', nrun=3, time.step=1, time.horizon=20, save.land=F, time.save=NA, out.path=NA )


ap.sbw <- function(scn, is.sbw = FALSE, is.harvesting = FALSE,
                       custom.params = NA, rcp = NA, 
                       nrun = 1, time.step = 1, time.horizon = 80, 
                       save.land = FALSE, time.save=5, out.path = NA, ...){
  
  
  ###########################################################################
  #################        0. LOAD INITIAL ARGUMENTS        #################
  
  #Library
  library(dplyr); library(stringr); library(raster); library(RANN); library(sp); library(raster); library(reshape); library(purrr)
  options(dplyr.summarise.inform=F)
  
  #Functions
  source("R/artificial.planting.r")
  source("R/buffer.mig.r")
  source("R/forest.transition.r")
  source("R/data.r")
  source("R/default.params.r")
  source("R/forest.mortality.r")
  source("R/harvest.area2.r")
  source("R/intens.def.curr.r")
  source("R/intensity.defoliation.r")
  source("R/neigh.influence.sbw.spread.r")
  source("R/sbw.outbreak2.r")
  source("R/spread.tonew.r")
  source("R/suitability.r")
  ## Function to select items not in a vector
  `%notin%` = Negate(`%in%`)
  
  #Data
  load("data/mask.rda")
  load("data/landscape.rda")
  load("data/prec_rcp85.rda")
  load("data/temp_rcp85.rda")
  load("data/prec_rcp45.rda")
  load("data/temp_rcp45.rda")
  load("data/default.tables.rda")
  
  
  #####################################################################
  #################        1. DATA PREPARATION        #################
  
  cat("A. Data preparation ...\n")
  
  ### 1.1. TIME STEPS ###
  ## Build the discrete time sequence according to time.step
  time.seq <- seq(0, time.horizon, time.step) 
  sbw.schedule <- c(0,35,70) #*** es pot eliminar?
  
  ## Set out.seq to save function outputs
  if(save.land){
      out.seq <- seq(0, time.horizon, time.save) 
      if(!all(out.seq %in% time.seq)){warning("Not all time steps in the output sequence provided are simulated.", call.=F)}
      if(is.na(out.path)) stop("Directory path to save outputs not provided") } 
  if(!(save.land)) {out.seq=0}

  
  ### 1.2. CELL RESOLUTION ###
  km2.pixel <- raster::res(mask)[1] * raster::res(mask)[2] / 10^6
  
  ### 1.3. PARAMETERS ###
  ## Get the list of default parameters and update user-initialized parameters
  params <- default.params()
  if(!is.na(custom.params)){
    # Check class of custom.params
    if((!inherits(customParams, "list"))) {
      stop("'custom.params' must be a named list")
    }
    ## Check that the names of the customized parameters are correct
    if(!all(names(custom.params) %in% names(params)))
      stop("Wrong custom parameters names")
    params <- custom.param
  }
  
  ### 1.4. SET CLIMATIC PROJECTIONS ###
  ## Load precipitation and temperature projections provided with the package according to the climatic scenario.
  if(!is.na(rcp)){
    prec.proj = get(paste0("prec_", rcp)); prec.chg = T
    temp.proj = get(paste0("temp_", rcp)); temp.chg = T
  }
  
  ### 1.5. INITIALIZE TRACKING DATA FRAMES ###
  breaks <- c(0,20,40,60,80,100,999)
  tags <- c("C10","C30", "C50", "C70", "C90", "OLD")
  track.spp.age.class <- data.frame(run=NA, year=NA, mgmt.unit=NA, spp=NA, age.class=NA, area=NA)
  track.suit.class <- data.frame(run=NA, year=NA, bioclim.domain=NA, potential.spp=NA, poor=NA, med=NA, good=NA)
  track.sbw.defol.intens = data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA)
  track.sbw.kill = data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA)
  track.cut <- data.frame(run=NA, year=NA, mgmt.unit=NA, spp=NA, age=NA)
  track.cut.sm <- data.frame(run=NA, year=NA, a.age=NA, ncells=NA)
  track.ap = data.frame(run=NA, year=NA, initial.spp=NA, ap.spp=NA, mgmt.unit=NA)
  track.ap.sm <- data.frame(run=NA, year=NA, ncells=NA)
  track.change = data.frame(cell.id=NA, run=NA, year=NA, spp=NA, age=NA, trans=NA)


  ######################################################################
  #################        2. START SIMULATIONS        #################
  
  cat("\n") 
  cat(paste("B. Simulations ...\n"))
  
  irun=1
  for(irun in 1:nrun){

    ## Main landscape data frame 
    tbls <- default.tables
    land <- landscape
    
    ## Upper truncate mean temperature of Black spruce, Yellow birch and Balsam fir
    land[land$spp == "EPN" & land$temp>= 2.4,]$temp <- 2.4
    land[land$spp == "BOJ" & land$temp>= 4.6,]$temp <- 4.6
    land[land$spp == "SAB" & land$temp>= 5.0,]$temp <- 5.0   
    
    ### SBW
    ## Decide (top-down) duration of the current outbreak, epidimic phase and collapse
    outbreak = 12 - params$current.duration  # testing
    collapse = params$collapse
    calm = params$calm
    preoutbreak = params$preoutbreak
    outbreak = rdunif(1,-1,1)+6+params$duration.last.outbreak - params$current.duration
    duration.last.outbreak = outbreak + params$current.duration
    phase = "outbreak"
    done = T
    
    ## Record initial distributions:
    ncell.def = c(NA, rep(sum(land$curr.intens.def>0),3))
    out = data.frame(run=irun, year=params$year.ini, phase=phase, group_by(land, curr.intens.def) %>% 
                       summarize(ncell=length(cell.id)) %>% mutate(pct=ncell/ncell.def))
    track.sbw.defol.intens = rbind(track.sbw.defol.intens, out)  

    
    # ## Determine the harvest regime for each species:  *** es pot eliminar?
    # ## even aged (1): 95% all conifers, PET, DB
    # ## uneven aged (0): 95% deciduous
    # ## other (2)
    # land <- mutate(land, rndm=runif(nrow(land)))
    # land$even <- 2
    # land$even[land$spp %in% c("EPN", "PET", "SAB", "OTH.RES.N", "OTH.RES.S", "OTH.FEU.N") & is.na(land$exclus) & land$rndm<=0.95] <- 1
    # land$even[land$spp %in% c("EPN", "PET", "SAB", "OTH.RES.N", "OTH.RES.S", "OTH.FEU.N") & is.na(land$exclus) & land$rndm>0.95] <- 0
    # land$even[land$spp %in% c("BOJ", "ERS", "OTH.FEU.S") & is.na(land$exclus) & land$rndm>0.95] <- 1
    # land$even[land$spp %in% c("BOJ", "ERS", "OTH.FEU.S") & is.na(land$exclus) & land$rndm<=0.95] <- 0    
    # 
    # ## Age classes distribution per species and management unit  *** per què aquí i no dins de t{}?? es pot eliminar?
    # land$age.class <- cut(land$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
    # track.spp.age.class <- rbind(track.spp.age.class, data.frame(run=irun, year=params$year.ini, 
    #                                                              filter(land, spp!="NonFor") %>% group_by(mgmt.unit, spp) %>% count(age.class) 
    #                                                              %>% mutate(area=n*km2.pixel)) %>% dplyr::select(-n))
    # ## Suitability classes distribution per bioclim.domain *** per què aquí i no dins de t{}??  es pot eliminar?
    # suitab <- suitability(land, params, tbls) 
    # aux <- left_join(suitab, dplyr::select(land, cell.id, bioclim.domain), by="cell.id") %>%
    #   group_by(bioclim.domain, potential.spp) %>% 
    #   summarise(poor=sum(suit.clim==0)*km2.pixel, med=sum(suit.clim==0.5)*km2.pixel, good=sum(suit.clim==1)*km2.pixel) 
    # track.suit.class <- rbind(track.suit.class, data.frame(run=irun, year=params$year.ini, aux))
    # rm(aux)
    
    
    ##############################
    ## Simulation per time step ##
    t=1
    for(t in time.seq){
      
      ## Print replicate and time step
      cat(paste0("Replicate ", irun, "/", nrun,". Time step: ", params$year.ini+t-time.step, "-", t+params$year.ini, "\n"))
      
      ###2.1. SET LAND TRANSITION AND LAND NEWSPP 
      land$transition=NA
      land$new_spp = NA
      
      ### 2.2. UPDATE CLIMATIC VARIABLES ###
      ## Update climatic variables at each time step if climate change is activated
      ## Column 1 is cell.index, the following columns are temp (precip) in 2020-2025, 2025-2030, 2030-2035, ... etc.
      ## The last column (temp95) then corresponds to the period 2095-2100
      ## The first time step (t=5) we start with climate 2020-2025
      if(temp.chg & t < time.horizon){
        cat("  Update temperature projections\n")
        aux <- temp.proj[,c(1,4+trunc(t/5))] 
        names(aux) <- c("cell.id", "temp")
        land <- dplyr::select(land, -temp) %>% left_join(aux, by="cell.id")
      }
      if(prec.chg & t < time.horizon){
        cat("  Update precipitation projections\n")
        aux <- prec.proj[,c(1,4+trunc(t/5))] 
        names(aux) <- c("cell.id", "prec")
        land <- dplyr::select(land, -prec) %>% left_join(aux, by="cell.id")
      }
      rm(aux)
      
      ### 2.3. UPDATE SUITABILITY
      suitab <- suitability(land, params, tbls)
      
      
      ##################################### PROCESSES OF CHANGE #####################################
      ### 1. SBW module

      if(is.sbw){
        cat ("  B.1. SBW \n")
        sbw.out <- integer()
        sbw.out <- sbw.outbreak(land,preoutbreak, outbreak, calm, collapse, duration.last.outbreak)
        
        kill.cells = sbw.out$kill.cells
        land.sbw = sbw.out$land.sbw    
        preoutbreak = sbw.out$preoutbreak 
        outbreak = sbw.out$outbreak    
        collapse = sbw.out$collapse    
        calm = sbw.out$calm 
        
        ## Tracking
        if(sum(land.sbw$curr.intens.def)>0){
          ncell.def = c(NA, rep(sum(land.sbw$curr.intens.def>0),length(unique(land.sbw$curr.intens.def))-1))
          track.sbw.defol.intens = rbind(track.sbw.defol.intens,
                     data.frame(run=irun, year=t+params$year.ini, phase=phase, group_by(land.sbw, curr.intens.def) %>%
                     summarize(ncell=length(cell.id)) %>% mutate(pct=ncell/ncell.def)))
        }
        else{
          track.sbw.defol.intens = rbind(track.sbw.defol.intens,
                     data.frame(run=irun, year=t+params$year.ini, phase=phase, curr.intens.def=0, ncell=nrow(land.sbw), pct=NA))
        }
        
        
        if(length(kill.cells)>0){
          aux = filter(land, cell.id %in% kill.cells) %>% group_by(spp, ny.def, curr.intens.def) %>% 
            summarise(area=length(spp)*km2.pixel) 
          names(aux) = c("spp", "ny.def", "curr.intens.def", "area")
          track.sbw.kill = rbind(track.sbw.kill, data.frame(run=irun, year=t+params$year.ini, aux))  
        }


        # Modifying land data frame
        land$age[land$cell.id %in% kill.cells] = 0
        land$tscomp = land$tscomp + params$time.step
        land$tssbw = land$tssbw + params$time.step   
        land$tssbw[land$cell.id %in% kill.cells] <- 0
        
        land$ny.def = land.sbw$ny.def
        land$ny.def[land$cell.id %in% kill.cells] = 0
        land$ny.def0 = land.sbw$ny.def0
        land$ny.def0[land$cell.id %in% kill.cells] = 0
        land$cum.intens.def = land.sbw$cum.intens.def
        land$cum.intens.def[land$cell.id %in% kill.cells] = 0
        land$curr.intens.def = land.sbw$curr.intens.def
        land$curr.intens.def[land$cell.id %in% kill.cells] = 0
        
      }#end is.sbw
      
       
      ### 2. Harvest area
      if(is.harvesting){
      cat ("  B.2. Harvest area \n")
      harv.out = integer()
      harv.out = harvest.area(land, params, default.tables, scn) 
      
      ## Tracking
      if(nrow(harv.out)>0)
        track.cut <- rbind(track.cut, data.frame(run=irun, year=t+params$year.ini, mgmt.unit=harv.out$mgmt.unit, spp=harv.out$spp, age=harv.out$age))
      if(nrow(harv.out)>0)
        track.cut.sm <- rbind(track.cut.sm, data.frame(run=irun, year=t+params$year.ini, a.age=trunc(mean(harv.out$age)), ncells=nrow(harv.out)))
      
      ## Modifying land data frame
      land$age[land$cell.id %in% harv.out$cell.id] <- 0
      land$age.class[land$cell.id %in% harv.out$cell.id] <- "C10" ##**** Class C10??
      land$new_spp[land$cell.id %in% harv.out$cell.id] <- NA
      land$transition[land$cell.id %in% harv.out$cell.id] <- "Harv"
      
      
      ### 3. Artificial planting
      cat ("  B.3. Artificial planting \n")
      ap.out = integer()
      ap.out <- artificial.planting(land, harv.out$cell.id, suitab, tbls, scn)
      
      ## Tracking
      if(nrow(ap.out)>0)
        track.ap <- rbind(track.ap, data.frame(run=irun, year=t+params$year.ini, initial.spp=ap.out$spp, ap.spp=ap.out$new.spp, mgmt.unit=ap.out$mgmt.unit))
      if(nrow(ap.out)>0)
        track.ap.sm <- rbind(track.ap.sm, data.frame(run=irun, year=t+params$year.ini, ncells=nrow(ap.out)))
      
      ## Modifying land data frame
      land$new_spp[land$cell.id %in% ap.out$cell.id] <- ap.out$new.spp
      land$transition[land$cell.id %in% ap.out$cell.id] <- "Harv_AP"
      
      ## Clear cutted cells without artificial planting
      no.ap.cells = filter(land, land$cell.id %in% harv.out$cell.id & !land$cell.id %in% ap.out$cell.id)
      
      } #end is.harvesting
      
      
      ### 4. Forest transition
      cat ("  B.4. Forest transition \n")
      
      ## After harvesting
      if(nrow(no.ap.cells)>0)
        land$new_spp[land$cell.id %in% no.ap.cells$cell.id] <- forest.transition(land, no.ap.cells$cell.id, suitab, params, default.tables, type.trans="C")
        land$transition[land$cell.id %in% no.ap.cells$cell.id] = "Harv_ForTrans"
      ## After SBW outbreak
      if(is.sbw){
        if(length(kill.cells)>0)
        land$new_spp[land$cell.id %in% kill.cells] <- forest.transition(land, kill.cells, suitab, params, tbls, type.trans="O")
        land$transition[land$cell.id %in% kill.cells] = "SBW_FotTrans"}
      ## Natural succession of tree spp at every 40 years starting at Tcomp = 70 **Que és tscomp??
      chg.comp.cells <- filter(land, (age-age.matu) %in% seq(40,400,40) & tscomp>=70) %>% dplyr::select(cell.id)
      if(length(unlist(chg.comp.cells))>0)
        land$new_spp[land$cell.id %in% unlist(chg.comp.cells)] <- forest.transition(land, unlist(chg.comp.cells), suitab, params, type.trans="S")  
        land$transition[land$cell.id %in% unlist(chg.comp.cells)] = "NatTrans"
      
        
      ### 5. UPDATE LAND.df
      land$spp = as.factor(ifelse(is.na(land$new_spp),as.character(land$spp), land$new_spp))
      land$age = land$age + params$time.step 
      
        
      ##################################### TRACKING AND SPATIAL OUTS #####################################
      cat ("  B. Traking spatial outputs \n")
      ## Age classes distribution per species and management unit **different than age.class???      
      land$ageClass <- cut(land$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
      track.spp.age.class <- rbind(track.spp.age.class, data.frame(run=irun, year=t+params$year.ini, 
                                                                   group_by(land, mgmt.unit, spp) %>% count(age.class) %>%  
                                                                     mutate(area=n*km2.pixel)) %>% dplyr::select(-n))
      ## Suitability classes distribution per bioclim.domain   
      suitab <- suitability(land, params, tbls) 
      aux <- left_join(suitab, dplyr::select(land, cell.id, bioclim.domain), by="cell.id") %>%
        group_by(bioclim.domain, potential.spp) %>% summarise(poor=sum(suit.clim==0)*km2.pixel, 
                                                              med=sum(suit.clim==0.5)*km2.pixel, good=sum(suit.clim==1)*km2.pixel) 
      track.suit.class <- rbind(track.suit.class, data.frame(run=irun, year=t+params$year.ini, aux))
      rm(suitab); rm(aux)
      
      
      ## Track change
      track.change <- rbind(track.change, data.frame(cell.id=land$cell.id, run=irun, year=t+params$year.ini, spp=land$spp, age=land$age, trans=land$transition))
      
      ## If required, save landscape data frame at the time steps required
      if(save.land & t %in% out.seq){
        if(!file.exists(out.path))
          dir.create(file.path(out.path), showWarnings = T) 
        saveRDS(land, file=paste0(out.path, "/landscape_", irun, "run_", t, "t.rds"))
      }
      
 
      ### Cleaning global environment 
      rm(sbw.out, land.sbw, kill.cells, ap.out, chg.comp.cells, harv.out, no.ap.cells)
      
      
  cat(paste0("       end of Time step: ", params$year.ini+t-time.step, "-", t+params$year.ini,"  \n"))
  
    }#end t
    
  cat(paste0("       end of Replicate=",irun, "/", nrun,"  \n \n"))
  }#end irun
  
  cat("\n", "C. Build outputs...\n")
  res <- list(SppByAgeClass = track.spp.age.class[-1,],
              SuitabilityClass = track.suit.class[-1,],
              LandChange = track.change[-1,],
              ArtificialPlanting = track.ap[-1,],
              ArtificialPlanting.sm = track.ap.sm[-1,],
              Cuts = track.cut[-1,],
              Cuts.sm = track.cut.sm[-1,],
              SBWDefol = track.sbw.defol.intens[-1,],
              SBWKill = track.sbw.kill[-1,])
  return(res)
  ## If required, save res
  if(save.land & t %in% out.seq){
    if(!file.exists(out.path))
      dir.create(file.path(out.path), showWarnings = T) 
    saveRDS(res, file=paste0(out.path, "/ap.sbw_results.rds"))
  }

}#end ap.sbw function

