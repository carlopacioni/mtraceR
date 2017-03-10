## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")
#  library(devtools)
#  install_github("carlopacioni/mtraceR")

## ---- eval=FALSE---------------------------------------------------------
#  
#  # Install devtools from CRAN if you haven't already
#  install.packages("devtools")
#  
#  library(devtools)
#  build_github_devtools()
#  
#  # Restart R before continuing
#  install.packages("./devtools.zip", repos=NULL)
#  
#  # Remove the package after installation
#  unlink("./devtools.zip")
#  
#  library(devtools)
#  install_github("carlopacioni/mtraceR")
#  

## ------------------------------------------------------------------------
library(mtraceR)
# locate the example files
data.path <- system.file("extdata", package="mtraceR")
# run mtrace
test_microsats <- mtrace(heating=TRUE, nchain=4, burn.in=0.4, dir.in=data.path, 
                         dir.out=NULL, trim=TRUE, thin=10, save2disk=FALSE, 
                         bayesallfile="bayesallfile_ms_3popsSt.gz")

## ---- echo=FALSE, out.width=450, fig.cap="__Figure 1.__ Overimposed traces of seven replicates (one colour for each replicate) of the analysis of a microsatellite locus and density of posterior distribution"----
knitr::include_graphics("Bad_trace2.tif")

## ---- echo=FALSE, out.width=500, fig.cap="__Figure 2.__ Gelman's statistics for locus 2 of the traces in Figure 1. "----
knitr::include_graphics("Bad_Gelm_Plot.tif")


## ------------------------------------------------------------------------
test_microsats[[3]][[2]]


## ------------------------------------------------------------------------
# for brevity, we subset the table
test_microsats[[1]][, 1:4]


## ------------------------------------------------------------------------
# for brevity, we subset the table
test_microsats[[2]][, c(3, 9)]


## ---- echo=FALSE, out.width=450, fig.cap="__Figure 3.__ Overimposed traces (one colour for each replicate) of the analysis of a microsatellite locus and density of posterior distribution"----
knitr::include_graphics("Good_trace2.tif")

## ---- echo=FALSE, out.width=500, fig.cap="__Figure 4.__ Gelman's statistics for locus 2 of the traces in Figure 3.  "----
knitr::include_graphics("Good_Gelm_Plot.tif")


## ------------------------------------------------------------------------
 mod.comp <- BF(lmLs=c(-4862.85, -4860.58, -4863.08, -4887.25))
 mod.comp

## ------------------------------------------------------------------------
bsp_plots <- BSP(dir.in=data.path, skylinefile=NULL, dir.out=NULL, all.loci=TRUE,
             overall=TRUE, params=1, gen=1, mu=1, save2disk=FALSE)

## ---- echo=FALSE, out.width=600------------------------------------------
knitr::include_graphics("BSP_each_locus.tif")


## ------------------------------------------------------------------------
bsp_plots[[3]]

## ---- echo=FALSE, out.width=400------------------------------------------
knitr::include_graphics("Multi_BSP.tif")


