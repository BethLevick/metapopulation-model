## source functions
source( "metapopulation-functions_15-01-22.r" )
##########################################################################################################
## Set up environment for analysing season by season ##
# get together the stuff we need from the file
# generate matrix of unique ids and their occupancy status each season
burrow <- unique(dat$Unique.Code)	#vector of all unique id numbers
x <- colnames(dat)	#take all the data column names
seasons <- c(x[9], x[11], x[13], x[15], x[17], x[21])	#take all the unique seasons

tmat <- matrix( ncol=length(seasons), nrow=length(burrow) )
rownames(tmat) <- burrow
colnames(tmat) <- seasons

for(k in 1:length(burrow)){	#for each unique burrow
tmp <- dat[dat$Unique.Code==burrow[k],]	#subset data frame to that burrow

vec <- c(tmp[,9], tmp[,11], tmp[,13], tmp[,15], tmp[,17], tmp[,21] )

#if( nrow(tmp) > 1 ){	# to check only 1 row for each unique id, unhash to re run
#print ( burrow[k] )
#}
#push vector of occupancy status of burrow in each season to corresponding row of matrix
tmat[k,] <- as.numeric(vec)

}
## will get warnings due to coerced data ##
# convert to a data frame and add some additional columns
#trans <- as.data.frame(tmat)
trans <- tmat
rm(seasons,burrow)	#remove to ensure no confusion when running models
##########################################################################################################
## convert occupancy statuses to just binary of occupied or not
for(k in 1:ncol(trans)){			#do this for each column in the trans data frame
trans[,k] <- binaryOcc( trans[,k] )
}
##########################################################################################################
## to generate difference matrix between burrows
## make reference data frame of burrow names, lat and long
#nams <- unique(dat$Unique.Code)
## find the associated latitude and longitude from the main data frame
#burrow <- data.frame( lat = sapply( nams, FUN=getCoord, dframe=dat, coordcol=7, namecol=6 ),
#	long = sapply( nams, FUN=getCoord, dframe=dat, coordcol=8, namecol=6 ) )
#rownames(burrow) <- nams

## generate distance matrix between all lat long pairs associated with a burrow
## saved as an .RDa loaded below - unhash to rerun
#dmat <- matrix( ncol=length(nams), nrow=length(nams), NA )
#rownames(dmat) <- nams
#colnames(dmat) <- nams
#for(k in 1:length(nams)){
#x <- burrow[nams[k],]
#	for(i in 1:length(nams)){
#	y <- burrow[nams[i],]
#	val <- dist( rbind(x,y) )
#	dmat[k,i] <- val
#	}
#}

#save( dmat, file="distMat.RDa" )
load( "distMat.RDa" )

### Convert to cartesian coordinates ###

## directly estimate mu from the data
losses <- data.frame( genshift=c("onetwo", "twothree", "threefour", "fourfive", "fivesix"), occupied=c(rep(NA,5)),
	becomelost=c(rep(NA,5)), fractionlost=c(rep(NA,5)) )
for(k in 1:ncol(trans)-1){
season1 <- trans[,k]
twoseasons <- as.character( paste( trans[,k], trans[,k+1], sep="" ) )
losses$occupied[k] <- length(which(season1==1))
losses$becomelost[k] <- length(which(twoseasons=="10"))
losses$fractionlost[k] <- losses$becomelost[k]/losses$occupied[k]
}

sumfractionlost <- sum(losses$becomelost)/sum(losses$occupied)
##########################################################################################################
## run the longterm-setup file to set up the predictive factors for the transitions data set
#source("longeterm-setup.r")
## strip the bits we need from the factor_df data frame
#properties <- data.frame( burrow=factor_df$burrow, lat=factor_df$lat, sand=factor_df$sand, allpl=factor_df$allpl, dunes=factor_df$binfactor3 )
## save this
#save( properties, file="burrow-properties.RDa" )
## and to reload
load( "burrow-properties.RDa" )



