sky.path<-"C:/Users/30373314/Documents/Woylie/Genetics/DNA analysis/allele scores & analysis/Historical samples Helen/ReAnalysis 2013/Migrate/NoDates/pan//Skyline//Cipres_Data"
sky.path<-"C:/Users/30373314/Documents/Woylie/Genetics/DNA analysis/mtDNA/Hist/ReAnalysis 2013/MIGRATE/FourPopsSt/skyline/Cipres_Data"
library(data.table)
library(ggplot2)
library(gridExtra)

dir.in=sky.path
skylinefile=NULL
dir.out=NULL
all.loci=TRUE
overall=TRUE
params=1
save2disk=TRUE

#' @param dir.in The local folder containing skylinefile files (default: NULL)
#' @param skylinefile The name of the skylinefile (default: NULL)
#' @param dir.out The local path to store the results. If NULL (default) then
#'  \code{dir.out=dir.in}
#' @param params A vector of parameter numbers to be plotted
#' @param all.loci Whether the parameter should be plotted for each locus
#' @param overall Whether the parameter should be plotted for the sum over all 
#'   loci (when \code{nloci>1})
#' @param save2disk Whether to save results to disk (default: TRUE)
#' @return See details for 
#'   additional outputs when \code{save2dis=TRUE}.
#' @import data.table
#' @import ggplot2
#' @import gridExtra
#' 
BSP <- function(dir.in=NULL, skylinefile=NULL, dir.out=NULL,  
                all.loci=TRUE, overall=TRUE, params=1, save2disk=TRUE){

  #----------------------------------------------------------------------------#
  # Helper funstions
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
    p <- dlocus[1:sel.min, ggplotDT(x=Age, y=Parameter.value) + 
                  geom_line() +
                  geom_ribbon(data=dlocus[1:sel.min, ],  alpha=0.2,
                                aes(x=Age, ymax=Upper, ymin=Lower)) +
                  theme_classic() + 
                  ylab(paste("Parameter", par)) + xlab("Time")]
    return(p)
  }
  
  plot.loci <- function(par, data, loci){
    dpar <- data[J(par), ]
    setkey(dpar, Locus)
    p.loci <- lapply(loci, plot.locus, dpar, par)
    plots <- do.call(arrangeGrob,  p.loci)
    return(plots)
  }
  #----------------------------------------------------------------------------#
  

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

rl <- readLines(paste0(dir.in, "/", "skylinefile"))
lstart <- grep(rl, pattern = "^1\t1\t")[1]
h <- c("Locus", "Parameter-number", "Bin", "Age", "Parameter-value", 
               "Parameter-Frequency", "Standard-deviation", "Counts-per-bin", 
               "Autocorrelation-per-bin")
d<-data.table((read.table(paste0(dir.in, "/", "skylinefile"), 
               skip=lstart - 1, header = F, col.names = h)))
nloci <- length(d[, unique(Locus)])
if(nloci > 1) nloci <- nloci - 1
setkey(d, Parameter.number)
d[, Standard.error := Standard.deviation / sqrt(Counts.per.bin)]
d[, Upper := Parameter.value + 1.96 * Standard.error]
d[, Lower := Parameter.value - 1.96 * Standard.error]  

lplots.byparams <- lapply(params, plot.loci, data=d, loci=1:nloci)


