#' Run migrate-n from within R
#' 
#' This is a convenience function to run migrate-n from within R
#' 
#' @param runOptions Character vector to pre-append to the command being issued. 
#'    The default assume that you are running migrate in parallel using openMPI
#'     (`"mpirun --use-hwthread-cpus"`)
#' @param migratePath Character vector with the path to migrate excutable
#' @param migrateExecutable The name of the executable (default `migrate-n-mpi`)
#' @param nProcessors The number of processor to be used. If `"auto"`, the maximum 
#'     number of processors available will be used
#' @param pathParmfile The path to the param file (default `"./"`
#' @param parmfile The name of the parmfile (default `"parmfile"`)
#' @export
runMigrate <- function(runOptions="mpirun --use-hwthread-cpus", 
                       migratePath="", migrateExecutable="migrate-n-mpi",
                       nProcessors="auto", 
                       pathParmfile="./", parmfile="parmfile") {
  oldwd <-getwd()
  on.exit(setwd(oldwd))
  setwd(pathParmfile)
  if(nProcessors == "auto") nProcessors <- parallel::detectCores()
  cmd <- paste(runOptions, "-np", nProcessors,
               file.path(migratePath, migrateExecutable), parmfile, "-nomenu")
  system(cmq)
}