## data set up functions
#########################################################################
###################### Data set up functions#############################
#########################################################################
## binary Occ ##
## Recode occupancy data as 1 and 0 from 0-4 categories ##
## function to convert a column of occupancy data into 1 and 0 where 1 is occupied and 0 is not
binaryOcc <- function(statusvec){
## take value in column away from 5, so will return 1 if occupied
new <- 5 - as.numeric(statusvec)
## where not equal to 1, replace with 0 (i.e. not occupied)
new[new!=1] <- 0
return(new)
}

## Example usage - best to use through sapply
# dat$bo.2011 <- binaryOcc( dat$occ.sp2011 )

#########################################################################
### getCoord ###
## function to get coordinates from name of location
## takes name of location, data frame to take from, column number of coordinate and column number of names
getCoord <- function( lname, dframe, coordcol, namecol ){
tmp <- dframe[dframe[,namecol]==lname,]
#print( tmp )
coord <- unique(tmp[,coordcol])[1]
return(coord)
}

#########################################################################
################ Functions to  estimate Parameters ######################
#########################################################################
#########################################################################

### findKi ###
## where occ.vec is the vector of occupancy statuses at the current generation
## dmat is the distance matrix
## params is a vector of values for C0 and L (as c(C0,L))
## and ref is the position of the burrow in the occupancy status vector
findKi <- function( occ.vec, dmat, params, ref ){
	## isolate the row relating to the burrow at position ref
	## and isolate from columns where the relating burrow is occupied
	## temporarily rewrite any NAs in the occupancy vector as 0
	occ.vec[is.na(occ.vec)] <- 0
	dists <- dmat[ref,occ.vec==1]
	#print(length(dists))
	## estimate and return Ki
	## where ki is the sum of C0 * ((exp(1)) ^( (-(abs(ri-rj)))/L ) )
	## for each pair of burrow ri compared to all other occupied burrows rj
	#print( -abs(dists) )
	#print( (exp(1)) ^( (-(abs(dists)))/params[2] ) )
	#print( params[1]*( (exp(1)) ^( (-(abs(dists)))/params[2] ) ) )
	return( sum( params[1]*( (exp(1)) ^( (-(abs(dists)))/params[2] ) ), na.rm=T ) )
}

## Example usage
# ki <- findKi( occ.vec=trans[1,], dmat, params=c(C0,L), ref=1 )

### ki2Ci ###
## convert Ki rates to Ci probabilities
ki2Ci <- function( ki ){
return( abs(ki)/(1+(abs(ki))) )
}

## example usage
# civec <- ki2Ci( kivec )

### nllKI ###
## where params is a vector of starting values for the parameters
## as c(C0,L) - find estimates using findKi
## dat is a matrix of empirical data to fit parameters to
## in the format of burrows as rows and time points as columns
## dmat is the distance matrix
nllKi <- function(params, dat, dmat){
## set nll to 0
nll <- 0
	## to estimate C0 and L
	for(i in 1:(ncol(dat)-1)){
	nll.gen <- 0
		for(k in 1:length(dat[,i])){
		## if there is a recording for the burrow
		if( !is.na(dat[k,i]) ){
		## if the burrow is current unoccupied, we want ki
			if( dat[k,i]==0 ){
			ki <- findKi( dat[,i], dmat, params, i )
			Ci <- ki2Ci( ki )
			#print(ki)
				## if there is a recording for the same burrow the next generation
				if( !is.na(dat[k,i+1]) ){
				## if in the next generation this burrow becomes occupied
				if( dat[k,i+1]==1 ){
				nll.gen <- nll.gen + log( Ci )
				}
				## or if it remains unoccupied
				if( dat[k,i+1]==0 ){
				nll.gen <- nll.gen + log(1- Ci)
				}
				}
			}
		}
		#print( nll.gen )
		}
	nll <- nll + nll.gen
	#print(nll)
	}
#print(-nll)
return(-nll)
}

## Example usage
## nll <- nllKi( params=c(C0,L), dat=trans, dmat )


### nllMu ###
## where params is a vector of length 1 of a starting value for the parameters
## dat is a matrix of empirical data to fit parameters to
## in the format of burrows as rows and time points as columns
nllMu <- function(params, dat){
# set nll to 0
nll <- 0
	## find nll given starting parameter at each generation
	for(i in 1:(ncol(dat)-1)){
	## set an empty variable to store the generations nll
	nll.gen <- 0
		## then for each record in this generation
		for(k in 1:length(dat[,i])){
			# if the data point is not missing
			if( !is.na(dat[k,i]) ){
				# if the burrow is occupied at this point in time
				# and the record is not missing the next generation
				if( dat[k,i]==1 && !is.na(dat[k,i+1]) ){
					# if the burrow becomes unoccupied
					if( dat[k,i+1]==0 ){
					nll.gen <- nll.gen + log( params[1] )
					#print( "+par" )
					}
					# if the burrow is occupied at this point in time
					# and still occupied at the next
					if( dat[k,i+1]==1 ){
					nll.gen <- nll.gen + log(1 - params[1] )
					#print( "+1-par" )
					}
					
				}
			
			}
		}
	## add the logged generational nll to the total
	nll <- nll + nll.gen
	}
## return the negative of the total
return( -(nll) )

}

## Example usage
# nll <- nllMu(params=c(0.1), dat=trans)



### generateParams
generateParams <- function(params, dat, dmat, send="Ki"){
	## if we request C0 and L using send="Ki" (default)
	if( send=="Ki" ){
	## estimate C0 and L, and then mu, and concatenate together as a vector
	opt <- optim( fn=nllKi, par=params, gr=NULL, dat, dmat, method="Nelder-Mead")#
	}#
	if(send=="Mu"){
	## otherwise (send=="Mu") estimate mu
	opt <- optim( fn=nllMu, par=params, gr=NULL, dat, method="Brent", lower=0, upper=1)
	}
return(opt)
}

## Example usage to estimate C0 and L (Where first in vector is C0 and second in vector is L)
## pars <- generateParams( params=c(0.01,0.01), dat=trans, dmat, send="Ki")
## Example usage to estimate mu
## mu <- generateParams( params=0.01, dat=trans, dmat, send="Mu")




#########################################################################	
nextGen <- function( prev, pars=c(mu,C0,L), dmat ){
nxt <- vector()
	## for each burrow in the vector
	for(i in 1:length(prev)){
	## if the burrow at this position is not missing
	if( !is.na(prev[i]) ){
		## if the burrow is occupied in the previous generation
		if( prev[i]==1 ){
		## it becomes unoccupied with the probability 1-mu
			x <- runif(1)
			#print(1 - (x + (pars[1])))
			#print( (x + (1-pars[1])) )
			nxt[i] <- floor( (x + (1-pars[1])) )
			#nxt[i] <- floor(x + (1-pars[1]))
			#print( floor(x + pars[1]) )
			#nxt[i] <- sample(x=c(1,0),size=1,prob=c(1-pars[1],pars[1]))
		}
		## if the burrow is unoccupied in the previous generation
		if( prev[i]==0 ){
		## it becomes occupied with the probability ki
			Ki <- findKi( prev, dmat, params=c(pars[2],pars[3]), ref=i)
			#print( paste("Ki:",Ki) )
			Ci <- ki2Ci( Ki )
			#print( paste("Ci:", Ci) )
			nxt[i] <- floor(runif(1) +  Ci)
			#nxt[i] <- sample(x=c(1,0),size=1,prob=c(ki2Ci( findKi( prev, dmat, params=c(pars[2],pars[3]), ref=i)),
				#1-ki2Ci( findKi( prev, dmat, params=c(pars[2],pars[3]), ref=i))) )
		}
	}else{
	## in some circumstances, burrows may be missing in the seed data
	## in this instance, put an empty burrow in this position in the next generation
		nxt[i] <- 0
		
	}
	}
return(nxt)
}

## Example usage
## x <- nextGen( seed, pars=c(mu,C0,L), dmat )	
	
#### Functions to run the simulation #######
metaSim <- function( seed, nruns, ngens, dmat, pars=c(mu,C0,L) ){
## set up a list to store results 
res <- list()
## for each run of the model
for(i in 1:nruns){
## set up a matrix to store results, with the seed as the first column
## col per generation, row per burrow
sim <- matrix( ncol=ngens, nrow=length(seed) )
sim[,1] <- seed
	## for each generation shift (gens-1)
	for(k in 2:ngens){
	## values in the column at k are determined by the column at position k-1
	## so vals in column k, where k-1==1 (occupied in previous generation)
	## become extinct with the probability mu
	sim[,k] <- nextGen( sim[,k-1], pars=pars, dmat )
	}
res[[i]] <- sim
}
return(res)
}

## Example usage
# results.list <- metaSim( seed, nruns, ngens, dmat, pars=c(mu,C0,L) )

######################################################################################
### plotModel ###
## function wrapper for commands to produce plot of model output
## just takes results list as argument
plotModel <- function( results.list, nruns, plotts=FALSE ){
results.summ <- matrix( ncol=ngens, nrow=nruns )
rownames(results.summ) <- seq(1,nruns,by=1)
colnames(results.summ) <- seq(1,ngens,by=1)

for(k in 1:nruns){
## for each run of the model (data frame position in the list)
## isolate the data frame
tmp <- results.list[[k]]
	## for each generation (column in the isolated data frame)
	for(i in 1:ngens){
	## find the proportion of burrows in this vector that are occupied over all burrows
	## send back to row k (run no) col i (gen no) of the summary matrix
	results.summ[k,i] <- length(which(tmp[,i]==1))
	}
}

## convert to proportions
burrows <- nrow(results.list[[1]])
results.summ <- results.summ/burrows

## convert the data set to proportions to plot against
data.summ <- vector()
## find number of burrows occupied in each generation
for(k in 1:ncol(trans)){
data.summ[k] <- length(which(trans[,k]==1))
}
## then convert to a proportion of all burrows
data.summ <- data.summ/burrows

## generate a plot title
plot.title <- paste("Proportion of burrows occupied per generation, over nrun=", nruns, sep="")
## plot data
plot( 0,0, main=plot.title, xlab="Time point", 
	ylab="Proportion Occupied", ylim=c(0,1), type="n", xlim=c(1,ngens) )

for(k in 1:nrow(results.summ)){
lines( c(1:ngens), results.summ[k,], col="lightpink" )
}

if( plotts==TRUE ){
lines( c(1:ncol(trans)), data.summ, col="black" )
}

legend( "topleft", legend=c("Simulation", "Empirical data"), col=c("lightpink", "black"), pch=c(16,16), bty="n" )
}