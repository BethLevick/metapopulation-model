### Bethany Levick, University of Liverpool 2014 ###
### To run metapopulation simulations of gerbil occupancy ###
##########################################################################################################
### to run metapopulation model ###
## set up environment, starting with opening working directory and then the data file
setwd( "C:/Users/Bethany/Dropbox/PhD/occupancy-model/main_R" )
dat <- read.table( "transitions.csv", sep=",", header=T, stringsAsFactors=F, fill=TRUE )
## source code to set up data/environment and source functions
## will return warnings - this is due to data coercion in generating transitions data frame and is expected
source( "setup.r" )
##########################################################################################################
### Set up for simulation ###
## set controls for model
## number of iterations of the model
nruns <- 100
## number of generations each iteration runs for
ngens <- 6
## starting (seed) vector of occupancies
seed <- trans[,1]

## estimate parameters
pars <- generateParams( params=c(0.01,0.01), dat=trans, dmat, send="Ki")
C0 <- abs(pars$par[1])
L <- abs(pars$par[2])
pars <- generateParams( params=0.01, dat=trans, dmat, send="Mu")
mu <- pars$par[1]

## save parameters to stop reloading every time
pars <- c(mu, C0, L)
save( pars, file="params.vec" )

### Load previously estimated parameters ####
load( "params.vec" )
mu <- pars[1]
C0 <- pars[2]
L <- pars[3]
##########################################################################################################
## run simulation
## list where each data frame is a full simulation, each col of which is a generation
## there is a data frame for the number of runs (boots) entered as nruns
results.list <- list()
results.list <- metaSim( seed, nruns, ngens, dmat, pars=c(mu,C0,L) )
## save results to file
save(results.list, file=paste( "results", Sys.Date(), sep="_" ))
runname <- paste( "results", Sys.Date(), sep="_" )
## reload most recent results produced
#load( paste( "results", Sys.Date(), sep="_" ) )
## This only works if loading the data on the same day!
## save the runname above and then can load using the runname (i.e. if loading after an overnight run)
## - need to add a random number to allow for multiple runs on the same day (prevent overwriting)
load( runname )
##########################################################################################################
## plot output
plotModel( results.list, nruns, plotts=TRUE )