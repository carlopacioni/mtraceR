#' Clean up a BayessAllFile
#' 
#' This function attempts to restore a corrupted file. See details to know more
#' about when this may happen, and what problems this function can fix.
#' 
#' When a run is interrupted, for example when it crashes or gets killed because 
#' it reached the limit set by a walltime, some lines (generally the last one) 
#' may be incomplete. This can cause issue when reading the file.
#' Another situation where there might be issues with the BayesAllFile is when 
#' you recover the run. Every time you restart a run, Migrate-n appends the headings 
#' and other text at the end of the file before starting logging again, which 
#' when read as a table, causes problems.
#' Also, when you restart the run multiple times, this text is appended multiple
#' times but Migrate-n uses the second instance to restart the run, so you restart
#' the run always from the same point rather than from the last processed locus.
#' This function tries to clean up the file by removing incomplete lines and the
#' unecessary text, so that it can be read as a table, and 
#' if you restart the run, it will resume form the last processed locus.
#' 
#'  @inheritParams mtrace
cleanBayesAll <- function(dir.in, bayesallfile="bayesallfile.gz") {
  need2write <- FALSE
  suppressWarnings(
  rlbay <- readLines(file.path(dir.in, bayesallfile))
  )
  grp <- grep("@@@@@@@@", rlbay)
  
  # check how long each line is (i.e. how many columns)
  rlbay_length <- sapply(strsplit(rlbay, split = "\t"), length)
  tb <- table(rlbay_length)
  # work out the correct number of columns
  nColCorrect <- as.integer(names(tb[which.max(tb)]))
  
  # Fix it if there is an error
  if(sum(rlbay_length[(grp[1] + 1):length(rlbay)] == nColCorrect) < 
     length((grp[1] + 1):length(rlbay))) { # the number of lines if it is correct
    need2write <- TRUE
    # exclude the initial text
    rlbay_temp <- rlbay[(grp[1] + 1):length(rlbay)]
    # keep the complete lines
    rlbay_fixed <- rlbay_temp[rlbay_length[(grp[1] + 1):length(rlbay)] == nColCorrect]
    # add the initial text
    rlbay <- c(rlbay[1:(grp[1] + 1)], rlbay_fixed)
    message("Problem(s) found with the bayesallfile, re-writing it...")
    file.rename(file.path(dir.in, bayesallfile), file.path(dir.in, sub("\\.gz", "_orig.gz", bayesallfile)))
    writeLines(rlbay, file.path(dir.in, sub(pattern = "\\.gz", "", bayesallfile)))
    
    R.utils::gzip(file.path(dir.in, sub(pattern = "\\.gz", "", bayesallfile)))
  }
  return(invisible(rlbay))
}
             
             