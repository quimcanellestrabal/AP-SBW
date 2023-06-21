#' Suitability
#'
#' Determines the environmental suitability of forest stands according to temperature, precipiation and soil type constraints
#'
#' @param land A \code{landscape} data frame with forest stand records in rows
#' @param params A list of default parameters generated by the function \code{default.params()} or a customized list of model parameters
#' @param tbls A list of default input tables as in \code{data(default.tables()} or a customized list of input tables
#' @return A data frame with enviornmental suitability ([0,1]) of forest stands
#' 
#' @export
#' 
#' @examples
#' data(landscape)
#' params = default.params()
#' suitability(landscape[runif(10,1,nrow(landscape)),], params, default.tables)
#' 

suitability <- function(land, params,tbls){
  
  # Vector with Potential Species
  potential.spp <- levels(land$spp)
  potential.spp <- c(potential.spp[str_length(potential.spp)==3], "OTH")
  
  # Compute soil and climatic suitability per SppGrp  
  # Final suitability corresponds to the minimum value between soil and climate suitability
  dta <- data.frame(cell.id=NA, potential.spp=NA, suit.soil=NA, suit.clim=NA)
  for(ispp in potential.spp){
    th.temp <- filter(tbls$temp.suitability, spp==ispp)[-1] 
    th.prec <- filter(tbls$prec.suitability, spp==ispp)[-1] 
    th.soil <- filter(tbls$soil.suitability, spp==ispp)[-1] 
    aux <- data.frame(cell.id=land$cell.id, potential.spp=ispp, temp=land$temp, prec=land$prec, soil=land$soil.type) %>%
           mutate(class.temp=as.numeric(cut(temp, th.temp)),
                  class.prec=as.numeric(cut(prec, th.prec)),
                  suit.temp=ifelse(is.na(class.temp), 0, ifelse(class.temp==2, 1, params$suboptimal)), #*segueixen sortint NA
                  suit.prec=ifelse(is.na(class.prec), 0, ifelse(class.prec==2, 1, params$suboptimal)),
                  suit.soil=as.numeric(th.soil[match(soil, c("T","O","R","S","A"))]),
                  suit.clim=pmin(suit.temp, suit.prec)) %>%
           dplyr::select(cell.id, potential.spp, suit.soil, suit.clim)
    dta <- rbind(dta, aux)
  }
  
  ## Upgrade climatic and soil suitability for "other" and "NonFor" SppGrp
  # subland$SuitClim[subland$PotSpp=="other"] <- 1
  dta <- rbind(dta, data.frame(cell.id=land$cell.id, potential.spp="NonFor", suit.soil=1, suit.clim=1))
  
  # Remove first NA row and order by cell.id
  dta <- dta[-1,]
  dta <- dta[order(dta$cell.id),]
  return(dta)
}