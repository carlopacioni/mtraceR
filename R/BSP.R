#------------------------------------------------------------------------------#
# Helper funstions
#------------------------------------------------------------------------------#

#' Helper function to plot ggplot directly from data.table.
#' 
#' See \url{from http://r.789695.n4.nabble.com/ggplot2-and-data-table-td3743483.html}
#'   for details
#'   
#' Credit to Timothée Carayol
#' @param ... Arguments to be passed (explicitly) to ggplot aesthetic
#' @export
aesDT <- function(...) {
  aes <- structure(list(...),  class = "uneval")
}

#' Helper function to plot ggplot directly from data.table.
#' 
#' See \url{from http://r.789695.n4.nabble.com/ggplot2-and-data-table-td3743483.html}
#'   for details
#'   
#' Credit to Timothée Carayol
#' @param ... ggplot elements
#' @export
ggplotDT <- function(...) {
  ggplot(, aesDT(...))
}

#------------------------------------------------------------------------------#


#' Generate Bayesian skyline plots
#' 
#' \code{BSP} generates Bayesian skyline plots for each analysed locus (when 
#'   relevant) and the sum over all loci. Plots are returned as ggplot objects 
#'   and, when \code{save2disk==TRUE}, pdf files. 
#'   
#' \code{BSP} will also return the data used to generate the BSPs as a data.table
#'   object and, when \code{save2disk==TRUE}, a csv file. These data are a 'clean
#'   up' version of these contained in the skyline file along with the upper and 
#'   lower limits. The column "Time" is the rescaled "Age" (from migrate output) 
#'   when \code{gen} and \code{mu} are provided. If you want to obtain the effective 
#'   population size (recall that Migrate-n returns theta=xmuNe for population size,
#'   where x is the number of ploid - e.g. for mtDNA, 4 for nuclear loci, etc. - 
#'   and mu is the mutation rate,
#'   if this is the parameter that you are considering) you could obtain this by 
#'   calculating Ne=theta/(x * mu). 
#'   
#' @param dir.in The local folder containing skylinefile files (default: NULL)
#' @param skylinefile The name of the skylinefile (default: NULL)
#' @param dir.out The local path to store the results. If NULL (default) then
#'  \code{dir.out=dir.in}
#' @param params A vector of parameter numbers to be plotted
#' @param gen The generation time in unit of time (e.g. years)
#' @param mu Mutation rate expressed either  'per generation' basis or in the 
#'  same unit as \code{gen} (e.g. years)
#' @param mu_unit The unit used for the mutation rate. Options are either 
#'    "generation" or "time"
#' @param all.loci Whether the parameter should be plotted for each locus
#' @param overall Whether the parameter should be plotted for the sum over all 
#'   loci (when \code{nloci>1})
#' @param save2disk Whether to save results to disk (default: TRUE)
#' @return A list where the first element is a data.table with skyline data and 
#'   the others are the BSPs.
#' @import data.table
#' @import ggplot2
#' @import gridExtra
#' @export
BSP <- function(dir.in=NULL, skylinefile=NULL, dir.out=NULL, all.loci=TRUE, 
                overall=TRUE, params=1, gen=1, mu=1, mu_unit="generation", 
                save2disk=TRUE) {

  #----------------------------------------------------------------------------#
  # Helper functions
  #----------------------------------------------------------------------------#
  #############################################
  # from http://r.789695.n4.nabble.com/ggplot2-and-data-table-td3743483.html
  # credit to Timothée Carayol     
  
  aesDT <- function(...) {
    aes <- structure(list(...),  class = "uneval")
  }
  
  ggplotDT <- function(...) {
    ggplot(, aesDT(...))
  }
  ##############################################
  
  plot.locus <- function(i, dpar, par){
    dlocus <- dpar[J(i), ]
    neg <- dlocus[, which(Standard.deviation < 0)[1]]
    infin <- dlocus[, which(Standard.deviation == Inf)[1]]
    sel.min <- suppressWarnings(min(neg, infin, na.rm = TRUE))
    if(sel.min == Inf) {sel.min <- dlocus[, length(Bin)]
        } else{
          sel.min <- sel.min - 1
        }
    p <- dlocus[1:sel.min, ggplotDT(x=Time, y=Parameter.value) + 
                  geom_line() +
                  geom_ribbon(data=dlocus[1:sel.min, ],  alpha=0.2,
                                aes(x=Time, ymax=Upper, ymin=Lower)) +
                  theme_classic() + 
                  ylab(paste("Parameter", par)) + xlab("Time") +
                  ggtitle(paste("Locus", i))]
    return(p)
  }
  
  plot.param <- function(par, data, loci){
    dpar <- data[J(par), ]
    setkey(dpar, Locus)
    p.loci <- lapply(loci, plot.locus, dpar, par)
    return(p.loci)
  }
  
  make.Grob <- function(lp.par){
    p.par <- do.call(arrangeGrob,  lp.par)
    return(p.par)
  }
  
  reduce <- function(lp.par){
    p <- lp.par[[1]]
    return(p)
  }
  
  plot2disk <- function(i, plots, params, suffix, dir.out) {
    ggsave(plots[[i]], dpi=600, 
           filename=paste0(dir.out, "/", "Parameter", params[i], suffix, ".pdf"))
    plot <- plots[[i]]
    save(plot, file=paste0(dir.out, "/", "Parameter", params[i], suffix, ".rda"))
  }
  #----------------------------------------------------------------------------#
  
  if(mu_unit == "time") mu <- mu * gen else
    if(mu_unit!="generation") stop("The argument 'mu_unit' can only be eithe 'generation' or 'time'")

  if(is.null(skylinefile) & is.null(dir.in)) {
    message("Please, select the skylinefile to import")
    full.path <- file.choose()
    skylinefile <- basename(full.path)
    dir.in <- dirname(full.path)
  } else {
    if(is.null(skylinefile)) skylinefile <- "skylinefile"
    if(is.null(dir.in)) 
      dir.in <- choose.dir("Select the folder where skylinefile is")
  }
  
  if(is.null(dir.out)) dir.out <- dir.in
  
  ################################################################
  # This works only with multilocus data, probably because of the end of line issue
  # d <- fread(paste0(dir.in, "/", "skylinefile"))
  # h <- c("Locus", "Parameter-number", "Bin", "Age", "Parameter-value", 
  #        "Parameter-Frequency", "Standard-deviation", "Counts-per-bin", 
  #        "Autocorrelation-per-bin")
  # setnames(d, h)
  ################################################################
  
  rl <- readLines(paste0(dir.in, "/", skylinefile))
  lstart <- grep(rl, pattern = "^1\t1\t")[1]
  h <- c("Locus", "Parameter-number", "Bin", "Age", "Parameter-value", 
                 "Parameter-Frequency", "Standard-deviation", "Counts-per-bin", 
                 "Autocorrelation-per-bin")
  d<-data.table((read.table(paste0(dir.in, "/", skylinefile), 
                 skip=lstart - 1, header = FALSE, col.names = h)))
  loci <- d[, unique(Locus)]
  all <- loci[which.max(loci)]
  if(length(loci) > 1) loci <- loci[-which.max(loci)]
  setkey(d, Parameter.number)
  d[, Upper := Parameter.value + 1.96 * Standard.deviation]
  d[, Lower := Parameter.value - 1.96 * Standard.deviation]  
  d[, Time := Age * gen / mu]
  if(save2disk == TRUE) write.csv(d, paste0(dir.out, "/", "sky.data.csv"),  
                                  row.names=FALSE)
  
  if(all.loci == TRUE) {
    lplots.eachLocus.byParams <- lapply(params, plot.param, data=d, loci=loci)
    lplots.byParams.Grob <- lapply(lplots.eachLocus.byParams, make.Grob)
    if(save2disk == TRUE) {
      i <- seq_along(params)
      lapply(i, plot2disk, plots=lplots.byParams.Grob, params, 
             suffix=".eachLocus", dir.out)
    } 
  }
  
  if(overall == TRUE) {
    lplots.overall.byParams <- lapply(params, plot.param, data=d, loci=all)
    lplots.overall.byParams <- lapply(lplots.overall.byParams, reduce)
    if(save2disk == TRUE) {
      i <-  seq_along(params)
      lapply(i, plot2disk, plots=lplots.overall.byParams, params, 
             suffix=".overall", dir.out)
    }
  }
  return(list(sky.data=d,
              if(all.loci == TRUE) lplots.byParams.Grob,
              if(overall == TRUE) lplots.overall.byParams))

}

#' Sum parameter values from different parameters
#' 
#' \code{meta.BSP} takes as input data equal to the first element of the output
#'   from \code{BSP} and combines (i.e. sum) the values of the parameters passed 
#'   with the argument \code{params}. This could be useful, for example, when 
#'   the analysis includes populations from a meta-population system and the user
#'   wants to calculate the size of meta-population (i.e. the total).
#'   
#' When multiple loci are used in the analysis, the argument \code{locus} can 
#'   be used to indicate which locus should be used. If \code{locus="max"},
#'   then the largest locus number is used (which corresponds to the sum over 
#'   all loci for multi-locus analyses). 
#'   
#' @param data The data generated by \code{BSP} 
#' @param dir.out The local path to store the results. If NULL (default) then
#'  \code{dir.out=wd()}
#' @param params A vector of parameter numbers to be combined and plotted
#' @param locus The locus to be considered
#' @param save2disk Whether to save results to disk (default: TRUE)
#' @return A list with the combined data (as data.table) and a BSP.
#' @import data.table
#' @import ggplot2
#' @export
   
meta.BSP <- function(data=NULL, dir.out=NULL, params=NULL, locus="max", 
                     save2disk=TRUE) {
  #----------------------------------------------------------------------------#
  # Helper funstions
  #----------------------------------------------------------------------------#
  equal.length <- function(par, dlocus, sel.min) {
    dpar <- dlocus[J(par), ]
    dpar.eq <- dpar[1:sel.min, ]
  }
  #----------------------------------------------------------------------------#
  
  if(is.null(data)) stop("Error: no data provided!")
  if(is.null(params)) stop("Please, provide parameter number(s)")
  if(is.null(locus)) stop("Please, indicate locus")
  if(is.null(dir.out)) dir.out <- wd()
  h <- c("Locus", "Parameter-number", "Bin", "Age", "Parameter-value", 
         "Parameter-Frequency", "Standard-deviation", "Counts-per-bin", 
         "Autocorrelation-per-bin", "Upper", "Lower", "Time")
  if(!identical(names(data), make.names(h))) stop("Error, data headings are not correct")
  data <- data.table(data)
  if(locus == "max") locus <- data[, max(unique(Locus))]
  setkey(data, Locus)
  dlocus <- data[J(locus), ]
  neg <- dlocus[, which(Standard.deviation < 0)[1], by="Parameter.number"]
  infin <- dlocus[, which(Standard.deviation == Inf)[1], by="Parameter.number"]
  sel.min <- suppressWarnings(min(neg[, V1], infin[, V1], na.rm = TRUE))
  bins <- dlocus[, length(Bin), by="Parameter.number"]
  if(sel.min == Inf) {sel.min <- min(bins[, V1])
  } else{
    sel.min <- sel.min - 1
  }
  setkey(dlocus, Parameter.number)
  ldpars.eq <- lapply(params, equal.length, dlocus, sel.min)
  dpars.eq <- rbindlist(ldpars.eq)
  meta <- dpars.eq[, lapply(.SD, sum), 
                 .SDcol=c("Parameter.value", "Upper", "Lower"), 
                 by="Time"]
  p <- meta[, ggplotDT(x=Time, y=Parameter.value) + 
                geom_line() +
                geom_ribbon(data=meta,  alpha=0.2,
                            aes(x=Time, ymax=Upper, ymin=Lower)) +
                theme_classic() + 
                ylab("") + xlab("Time")]
  if(save2disk == TRUE) {
    pars <- paste0(params, collapse="-")
    write.csv(meta, paste0(dir.out, "/", "meta.", pars, ".csv"), row.names=FALSE)
    ggsave(p, dpi=600, filename=paste0(dir.out, "/", "meta.", pars, ".pdf"))
    save(p, file=paste0(dir.out, "/", "meta.", pars, ".rda"))
  }
  return(list(data=meta, plot=p))
}

#' Generate multiple Bayesian skyline plots within one figure
#' 
#' \code{multi.BSP} takes as 2 separate input data equal to the first element of 
#'   the output from \code{BSP} reduced to a unique parameter and a unique locus 
#'   and generates a figure where two BSP are concurrently plotted.
#'     
#' @param data1,data2 The data generated by \code{BSP} to be used for the BSPs
#' @param col1,col2 The colour to be used for the BSPs
#' @param y.title The title for the y-axis (default="")
#' @param dir.out The local path to store the results. If NULL (default) then
#'  \code{dir.out=getwd()}
#' @param save2disk Whether to save results to disk (default: TRUE)
#' @return A list with the combined data (as data.table) and a BSP.
#' @import data.table
#' @import ggplot2
#' @export

multi.BSP <- function(data1=NULL, data2=NULL, col1="black", col2="red",
                      y.title="", dir.out=NULL, save2disk=TRUE) {
  #----------------------------------------------------------------------------#
  # Helper funstions
  #----------------------------------------------------------------------------#
  err.handling <- function(data.x) {
    if(is.null(data.x)) stop("Error: no data provided!")
    h <- c("Parameter.value", "Upper", "Lower", "Time")
    if(sum(h %in% names(data.x)) < 4) stop("Error, headings are not correct")
    data <- data.table(data.x)
    return(data)
  }
  
  adj.length <- function(data.x){
    if("Standard.deviation" %in% names(data.x)) {
      neg <- data.x[, which(Standard.deviation < 0)[1]]
      infin <- data.x[, which(Standard.deviation == Inf)[1]]
      sel.min <- suppressWarnings(min(neg, infin, na.rm = TRUE))
      if(sel.min == Inf) {sel.min <- data.x[, length(Time)]
      } else{
        sel.min <- sel.min - 1
      }
      return(data.x[1:sel.min, ])      
    } else {
      return(data.x)
    }
    
  }
  
  write2disk <- function(i, data) {
    write.csv(data[[i]], paste0(dir.out, "/", "data.", i, ".csv"), row.names=FALSE)
  }
  
  #----------------------------------------------------------------------------#
        
  if(is.null(dir.out)) dir.out <- getwd()
  data <- list(data1, data2)
  data <- lapply(data, err.handling)
  data <- lapply(data, adj.length)
  
  p <-  ggplot(data[[1]], aes(x=Time, y=Parameter.value)) + 
                   geom_line(colour=col1) +
                   geom_ribbon(data=data[[1]], alpha=0.2, fill=col1,
                               aes(x=Time, ymax=Upper, ymin=Lower)) +
                   theme_classic() + ylab(y.title) + xlab("Time")
  
  pfin <- p + geom_line(data=data[[2]], aes(x=Time, y=Parameter.value), colour=col2) +
                      geom_ribbon(data=data[[2]], alpha=0.2, fill=col2,
                                  aes(x=Time, ymax=Upper, ymin=Lower)) +
                      theme_classic() + ylab(y.title) + xlab("Time")
    
    if(save2disk == TRUE) {
      i <- 1:2
      lapply(i, write2disk, data=data)
    ggsave(pfin, dpi=600, filename=paste0(dir.out, "/", "multi.BSP.pdf"))
    save(pfin, file=paste0(dir.out, "/", "multi.BSP.rda"))
  }
  return(list(data=data, plot=pfin))
}