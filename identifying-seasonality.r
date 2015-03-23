## Bethany Levick 2015
## To identify seasonal patterns in colonisation and extinction events
## set up environment, starting with opening working directory and then the data file
setwd( "C:/Users/Bethany/Dropbox/PhD/occupancy-model/main_R" )
dat <- read.table( "transitions.csv", sep=",", header=T, stringsAsFactors=F, fill=TRUE )
## source code to set up data/environment and source functions
## will return warnings - this is due to data coercion in generating transitions data frame and is expected
source( "setup.r" )

## dat is the original data set
## trans is the cleaned transitions data set
## first - do we get any new colonisation events occuring from an autumn to the following spring
## only have one instance of this, Au2011 to Sp2012
## paste them together to provide a two character description of occupancy events
ausp <- paste(trans[,3],trans[,4], sep="")

## find events where an unoccupied burrow becomes occupied
length(which(ausp=="01"))
> 25
## so there are colonisation events from Autumn to Spring

####################################################################
## look at the differences in mu by seasons
spsu <- paste(trans[,1], trans[,2], sep="" )
suau <- paste(trans[,2], trans[,1], sep="" )

## losses from spring to summer
length(which(spsu=="10"))
> 77
## losses from summer to autumn
length(which(suau=="10"))
> 52
## losses from autumn to spring
length(which(ausp=="10"))
> 99