#------------------------------------------------------------------------------#
# Package info
#------------------------------------------------------------------------------#
#' mtraceR: an R package to analyse migration-n output
#'
#' \code{mtraceR} can be used to visualise trace of MCMC chains, calculate 
#' Effective Sample Size (ESS), Gelman and Rubin's convergence diagnostic (and 
#' plots) or Bayes factor, and generate skyline plots.
#'
#' The main function \code{mtrace} is essentially a wrapper for the R package 
#' \code{coda} optimised for migrate-n outputs. It is intended to facilitate and 
#' speed up post-analysis data processing.
#' 
#' COMPLETE DESCRIPTION/DETAILS
#'
#' @section Documentations:
#' Use \code{help(package = "mtraceR")} for a list of \code{mtraceR} functions
#' and their specific documentations.
#'
#' A more detailed description of the package and functions can be opened with:
#' \code{vignette(package="mtraceR", topic="mtraceR-package")}.
#'
#' More vignettes may be come available in the future. Use
#' \code{vignette(package="mtraceR")} to see all the available vignettes.
#'
#' @section Citation:
#' If you use \code{mtraceR}, please cite:
#' ...
#' 
#'
#' @section Get in touch:
#' Please, use \url{https://github.org/carlopacioni/mtraceR/issues} to report
#' any issues with \code{mtraceR}. If unsure, or for feedback, contact us at:
#' carlo.pacioni 'at' gmail.com or anders???.
#'
#' @docType package
#'
#' @name mtraceR
NULL


#' Diagnostics of migrate-n MCMC chains
#' 
#' \code{mtrace} parses migrate-n output file (bayesallfile) and calculates the
#' Effective Sample Size (ESS) for each locus. If \code{save2dis=TRUE} these 
#' are also saved to disk together with a PDF with the trace and density plots 
#' for each locus and each replicate.
#' 
#' If there is more than 1 replicate, \code{mtrace} returns also Gelman and 
#' Rubin's convergence diagnostic for each locus. When \code{save2dis=TRUE} 
#' relevant diagnostic plots are also saved as a separate PDF.
#' 
#' When no bayesallfile name and directory are provided (default), the output 
#' file is selected with an iteractive window. Alternatively, when the directory 
#' path is provided with \code{dir.in} and \code{bayesallfile=NULL}, the latter 
#' is automatically set to "bayesallfile.gz". If the name of the file is provided,
#' but \code{dir.in=NULL}, the directory path is selected iteractively.
#' 
#' If \code{dir.out=NULL} (default) then \code{dir.out=dir.in}. 
#' 
#' @param heating Whether heating is used: (\code{TRUE} or \code{FALSE})
#' @param nchain Number of chains if \code{heating=TRUE}
#' @param burn.in Length of burn-in as proportion of MCMC (default: 0.1)
#' @param npop The number of populations in the analysis
#' @param estim.npar The number of estimated parameters (when a custom migration 
#' matrix is used)
#' @param trim Whether trim data to estimated parameters (default:  TRUE)
#' @param bayesallfile The name of the bayesallfile (default: NULL)
#' @param dir.in The local folder containing migrate-n output files (default: NULL) 
#' @param dir.out The local path to store the results. If NULL (default) then
#'  \code{dir.out=dir.in}
#' @param save2disk Whether to save results to disk (default: TRUE)
#' @return a list with the ESS for each locus, a summary of ESS across all loci 
#'   and Gelamn's diagnostic if more replicates were run. See details for 
#'   additional outputs when \code{save2dis=TRUE}.
#' @import data.table
#' @import coda
#' @export


mtrace <- function(heating=TRUE, nchain=4, burn.in=0.1, npop=1, estim.npar=NULL, 
                   trim=TRUE, bayesallfile=NULL, dir.in=NULL, dir.out=NULL, 
                       save2disk=TRUE) {
    
  #----------------------------------------------------------------------------#
  # Helper funstions
  #----------------------------------------------------------------------------#
  
     # Split data in bayesallfile by replicates
     # x Element of a list passed with lapply()
        byRepl <- function(x) {
          split(x, x$Replicate)
        }
       
     # Remove first two columns (Steps & Locus)
     # x Element of a list passed with lapply()
     noLocus <- function(x) {
       ncol <- dim(x[[1]])[2]
       noloc <- lapply(x, function(slim) slim[, 3:ncol])
       return(noloc)
     }
     
     # Convert object in mcmc
     # x Element of a list passed with lapply()
     make.mcmc <- function(x) {
       x.mcmc <- lapply(x, mcmc)
       return(x.mcmc)
     }
     
     # Remove burn-in
     # x Element of a list passed with lapply()
     # burn burn-in. This is calculate automatically based on burn.in
     rm.burn <- function(x, burn) {
       x.burn <- lapply(x, window, start=round(burn))
       return(x.burn)
     }
     
     # Trim down data by removing unnecessary columns
     # x Element of a list passed with lapply()
     trim.data <- function(x,  last.col) {
       trimmed <- lapply(x, function(tr) tr[, 7:(last.col - 3)])
       trimmed <- mcmc.list(trimmed)
       return(trimmed)
     }
     
     # Gelman and Rubin's convergence diagnostic
     # x Element of a list passed with lapply()
     diagnostic <- function(x) {
       x.diagn <- lapply(x, gelman.diag, autoburnin=FALSE)
       return(x.diagn)
     }
     
     # Save to disk Gelman and Rubin's convergence diagnostic
     # x A list where each element is the diagnostic for a locus
     # i A sequence integer vector of length equal to the number of loci 
     save.gel.diag <- function(i, x) {
       Locus <- rep(i, dim(x[[i]]$psrf)[1])
       Parameter <- rownames(x[[i]]$psrf)
       x.locus <- as.data.frame(x[[i]]$psrf)
       x.locus <- cbind(Parameter, Locus, x.locus)
       return(x.locus)
     }
     
     # Trace and density plots with locus name
     # x A list where each element is a mcmc.list with replicates for one 
     #  locus
     # i A sequence integer vector of length equal to the number of loci 
     plot.with.name <- function(i, x) {
       plot(x[[i]], ylab=paste("Locus", names(x)[i]), ask=FALSE)
     }
     
     # Plot of Gelman and Rubin's convergence diagnostic
     # x A list where each element is a mcmc.list with replicates for one 
     #  locus
     # i A sequence integer vector of length equal to the number of loci 
     gelman.plot.with.name <- function(i, x) {
       gelman.plot(x[[i]], ylab=paste("Locus", names(x)[i]), autoburnin=FALSE)
     }
     
  
  ###---------------------------- Data handling -----------------------------###
  # data.table 1.9.6 should work, 
  # data.table 1.9.5  though has an issues with the end of line 
  # so, for now, stick to read.table
  
  if(is.null(bayesallfile) & is.null(dir.in)) {
    message("Please, select the bayesallfile to import")
    full.path <- file.choose()
    bayesallfile <- basename(full.path)
    dir.in <- dirname(full.path)
  } else {
    if(is.null(bayesallfile)) bayesallfile <- "bayesallfile.gz"
    if(is.null(dir.in)) dir.in <- choose.dir("Select the folder where migate-n output files are")
  }
  
  if(is.null(dir.out)) dir.out <- dir.in
  
  message(paste("Parsing", bayesallfile))
  rl <- readLines(paste0(dir.in, "/", bayesallfile), n=100)
  patt <- grep(pattern = "^# @@@@@@@@", rl)
  
  h <- make.names(read.table(paste0(dir.out, "/", bayesallfile), header=F, skip=patt, nrow=1, 
                             colClasses="character"))
  data <- read.table(paste0(dir.out, "/", bayesallfile), header=F, skip=patt + 1, colClasses="numeric", 
                     comment.char="")
  tot.npar <- npop^2
  if(is.null(estim.npar)) estim.npar <- tot.npar
  cutoff <- 2 + tot.npar - estim.npar + if(heating == TRUE) nchain
  parlast <- length(h) - cutoff 
  h <- c(h, rep("EMPTY", tot.npar - estim.npar))
  names(data) <- h
  repl <- length(unique(data$Replicate))
  
  l.locus <- split(data[, 2:parlast], data$Locus)
  l.locus.repl <- lapply(l.locus, byRepl)
  burn <- burn.in * dim(l.locus.repl[[1]][[1]])[1]
  l.locus.repl <- lapply(l.locus.repl, noLocus)
  l.locus.mrepl <- lapply(l.locus.repl, make.mcmc)
  l.locus.mreplno.burn <- lapply(l.locus.mrepl, rm.burn, burn=burn)
  l.mlLocus.mrepl <- lapply(l.locus.mreplno.burn, mcmc.list)
  if(repl > 1 | trim == TRUE) {
    l.mlLocus.mrepltr <- lapply(l.mlLocus.mrepl, trim.data, last.col=parlast)
  }
  message("Done!")
  
  
  ###-------------------------- Gelman diagnostic ---------------------------###
  if(repl > 1) {
    g.diag <- lapply(l.mlLocus.mrepl, gelman.diag, autoburnin=FALSE, 
                     multivariate = F)
  }
  
  
  ###-------------------------------- ESS -----------------------------------###
  message("Calculating ESS after removal of burn-in...")
  ESSlist <- lapply(l.mlLocus.mrepl, effectiveSize)
  ESS <- as.data.frame(ESSlist)
  names(ESS) <- paste("Locus", 1:length(ESSlist))
  EES.summary <- apply(ESS, 1, summary)
  message("Done!")
  
  
  ###---------------------------- save2disk ---------------------------------###
  if(save2disk == TRUE) {
    message("Saving reults to disk...")
    write.csv(ESS, file=paste0(dir.out, "/", "ESS.csv"))
    write.csv(EES.summary, file=paste0(dir.out, "/", "EES.summary.csv" ))
    
    i <- seq_along(l.mlLocus.mrepl)
    
    pdf(file=paste0(dir.out, "/", "TracePlots_burninRemoved.pdf"))
    if(trim == TRUE) {
      lapply(i, plot.with.name, x=l.mlLocus.mrepltr)
    } else {
      lapply(i, plot.with.name, x=l.mlLocus.mrepl)
    }
    
    dev.off()
    
    if(repl > 1) {
      gelm.res <- lapply(i, save.gel.diag, x=g.diag)
      write.csv(data.table::rbindlist(gelm.res), 
                file=paste0(paste0(dir.out, "/", "gelm.diag.csv")))
      pdf(file=paste0(dir.out, "/", "gelm.plots.pdf"))
      lapply(i, gelman.plot.with.name, x=l.mlLocus.mrepltr)
      dev.off()
    }
    message("Done!")
  }
  
  ###------------------------------------------------------------------------###
  return(list(ESS=ESS, 
              EES.summary=EES.summary, 
              gelm.diag=if(repl > 1) g.diag else NULL))
}


BF <- function(lmLs=NULL) {
  # error handling
  if(!is.numeric(lmLs)) stop("Please, pass a numeric vector of log 
                             likelihood with the argument lmLs")
  ref.mod <- max(lmLs)
  LBF <- lmLs - ref.mod
  BF <- exp(LBF )
  mod.prob <- round(BF / sum(BF), digits=3)
  mod.rank <- rank( - LBF)
  return(data.frame(lmL=lmLs, LBF, mod.rank, mod.prob))
}