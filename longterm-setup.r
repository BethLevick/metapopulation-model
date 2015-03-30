## Bethany Levick, University of Liverpool 2014 ##
## Code to set up environment for running long term occupancy models from transitions file ##
################################################################################################################################
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
trans <- as.data.frame(tmat)

rm(seasons,burrow)	#remove to ensure no confusion when running models

################################################################################################################################
## generate new column of outcome variable - 
trans$sea_occ <- c(rep(0,(nrow(trans))))	#percentage of seasons occupied
trans$occ_class <- c(rep(0,(nrow(trans))))	#category of number of seasons occupied

## 11-03-13
## to use binomial error function, need to know number of times data recorded
trans$wgt <- c(rep(0,nrow(trans)))	#weights column
for(k in 1:nrow(trans)){			#for each row of data, take the string of numbers generated from each of the season columns
trans$wgt[k] <- length(which(trans[k,1:6]>-1))	#find length of vector where all in this row not NA
}

classes <- c("a","b","c","d","e","f","g")	#for categories, with a being 0 occupied and g being 6 occupied

for(k in 1:nrow(trans)){
	ovec <- c(trans[k,1:6])	#generate vector of occupancy statuses across 6 seasons
	x <- which( ovec==4 )	#find how many of these are occupied
	
	trans$occ_class[k] <- classes[(length(x)+1)]	#the class is in the position length+1, as length 0 is position 1
	
	po <- (length(x)/trans$wgt[k])	#find proportion of all these seasons that this burrow was occupied
	trans$sea_occ[k] <- po	#pass this to the relevant column in the data frame
	}
##################################################################################################################################
## set up new data frame with predictive factors from main data frame
factor_df <- data.frame( landscape=as.factor(dat$Bscharact_all) )
#copy across seasons occupied and weights
factor_df$sea_occ <- trans$sea_occ
factor_df$wgt <- trans$wgt

# copy across latitude
factor_df$lat <- dat$Latitude..N.

# copy across tree and format
factor_df$tree <- as.factor(dat$tree)			#presence or absence of tree
factor_df$tree[factor_df$tree=="no data"] <- NA	#replace "no data" with NA

# add burrows & sectors for random effect
factor_df$burrow <- as.factor(dat$Unique.Code)
factor_df$sector <- as.factor(dat$Sector)

#function to perform regex search on a character string of the regex reg
#where l is the string to search and reg is the exact regex e.g. "1" or "[56]"
#will return a 1 if present and a 0 if absent in the character string
#neg is the value to print if the regex is not found in the string
#pos is the value to print if the regex is found in the string
classByRegex <- function(cs, reg, neg, pos){
val <- neg										#initally set the value to as specified in neg
x <- grep( reg, cs, value=F, perl=T )			#perform regex search in string cs
if(length(x)>0){								#if regex found, then length of result will be >0
val <- pos										#in which case specify value to be as specified in pos
}
return(val)										#in either case, return the value
}

## example usage, for sand (class =1)
##factor_df$sand <- sapply( factor_df$landscape, FUN=classByRegex, reg="1", neg=0, pos=1 )

## in place of factor 1 (ordinal sand) and factor 4 (alluvial plain)
## separate descriptors for sand, loam and alluvial plain
factor_df$factor2 <- sapply( factor_df$landscape, FUN=classByRegex, reg="[56]", neg=0, pos=1 )
factor_df$binfactor3 <- sapply( factor_df$landscape, FUN=classByRegex, reg="[DEF87]", neg=0, pos=1 )
factor_df$sand <- sapply( factor_df$landscape, FUN=classByRegex, reg="1", neg=0, pos=1 )
factor_df$sandloam <- sapply( factor_df$landscape, FUN=classByRegex, reg="2", neg=0, pos=1 )
factor_df$clayloam <- sapply( factor_df$landscape, FUN=classByRegex, reg="3", neg=0, pos=1 )
factor_df$allpl <- sapply( factor_df$landscape, FUN=classByRegex, reg="4", neg=0, pos=1 )
factor_df$loam <- sapply( factor_df$landscape, FUN=classByRegex, reg="[23]", neg=0, pos=1 )

# recode season
factor_df$binseason <- sapply( factor_df$landscape, FUN=classByRegex, reg="SP2013", neg="not spring", pos="spring" )

# convert latitude to region, with 45.5 degrees north being the cut off for the north and south regions
factor_df$reg <- factor_df$lat
factor_df$reg[factor_df$reg>45.4] <- "N"
factor_df$reg[factor_df$reg!="N"] <- "S"

## ensure set as factor
factor_df$sand <- as.factor(factor_df$sand)
factor_df$sandloam <- as.factor(factor_df$sandloam)
factor_df$clayloam <- as.factor(factor_df$clayloam)
factor_df$allpl <- as.factor(factor_df$allpl)
factor_df$loam <- as.factor(factor_df$loam)
factor_df$factor2 <- as.factor(factor_df$factor2)
factor_df$binfactor3 <- as.factor(factor_df$binfactor3)
factor_df$binseason <- as.factor(factor_df$binseason)
factor_df$reg <- as.factor(factor_df$reg)

#############################################################################################################
## convert response to integer to allow glmm to run
#factor_df$sea_occ <- as.integer( factor_df$sea_occ )