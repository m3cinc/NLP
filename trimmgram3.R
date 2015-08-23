# trimmgram3.R
options(mc.cores=4)
memory.limit(size=28000)
setwd("V:/") 
Sys.info()[1:5]                         # system info but exclude login and user info
userdir<-getwd()
# user-defined startup Capstone directory
library(dplyr)
library(magrittr)
library(stringi)
library(pbapply)
library(compiler)
local <- TRUE                           # to avoid data redownload
datadir <- "./data"
datasetsdir <- paste0(datadir,"/datasets")
cleandir <- paste0(datadir,"/clean")
gramdir <- paste0(datadir,"/gram")
xgramdir <- paste0(datadir,"/xgram")
lgramdir <- paste0(datadir,"/lgram")
mgramdir <- paste0(datadir,"/mgram")
sgramdir <- paste0(datadir,"/sgram")
tgramdir <- paste0(datadir,"/tgram")
ngramdir <- paste0(datadir,"/ngram")
maxlen <- 6L
if (!file.exists(datadir)) { dir.create(datadir) }  #  data will reside in datadir
if (!file.exists(datasetsdir)) { dir.create(datasetsdir) }  #  datasets will reside in datasetsdir
if (!file.exists(gramdir)) { dir.create(gramdir) }  #  gram will reside in datadir
if (!file.exists(xgramdir)) { dir.create(xgramdir) }  #  xgram will reside in datadir
if (!file.exists(lgramdir)) { dir.create(lgramdir) }  #  lgram will reside in datadir
if (!file.exists(mgramdir)) { dir.create(mgramdir) }  #  mgram will reside in datadir
if (!file.exists(sgramdir)) { dir.create(sgramdir) }  #  sgram will reside in datadir
if (!file.exists(tgramdir)) { dir.create(tgramdir) }  #  tgram will reside in datadir
if (!file.exists(ngramdir)) { dir.create(ngramdir) }  #  ngram will reside in datadir

initialize <- cmpfun( function(x, z, n = maxlen, level = 0L) {
        #' load('CS50GramX0.RData')
        #' @param x the name of the grambase (e.g. "CS50")
        #' @param z the size of the grambase (eg) gramdir, xgramdir, lgramdir, mgramdir, sgramdir, tgramdir, ngramdir
        #' @param n integer the number of grams to load
        #' @param level integer to trim
        #' @return a list of grams in x 
        #' @examples 
        #' initialize(CS50Gram)   initializes up to 7-gram 
        #' initialize(CS50Gram,3) initializes up to trigram
        stopifnot (is.integer(n), is.finite(n), n > 0L, n <= maxlen,
                   is.integer(level), is.finite(level), level > -1L, is.character(x))
        extract <- magrittr::extract
        filelist <- list()
        filelist <- sapply(1:n, function(i) {filelist[i] <- paste0(x, "Gram", i, "0")})
        x <- list()
        print ("Reading Gram:")
        x <- sapply (seq_along (filelist), function(i) { 
                paste( z, filelist[i], sep = "/" ) %>% paste0(".RData") -> filename
                load(file = filename)
                print (c(filelist[i],", "))
                x[[i]] <- Vocabulary  })  # only retain the Gram
        if (level > 0) { x <- pblapply(1:n, function(i) { x[[i]] %<>% extract(.> level) } )
                         x <- lapply(1:n, function(i) { x[[i]] %<>% extract %>% -level } ) }
        x
}, options = list( optimize = 3) )

initializeX <- cmpfun(function(x, z, n = maxlen ,level = 0L) {
        #' load('CS50GramX0.RData')
        #' @param x the name of the grambase (e.g. "CS80") 
        #' @param z the size of the grambase (eg) gramdir, xgramdir, lgramdir, mgramdir, sgramdir, tgramdir, ngramdir
        #' @param n integer the number of grams to load
        #' @param level integer to trim
        #' @return a list of grams in x 
        #' @examples 
        #' initialize(CS50Gram)   initializes up to 7-gram 
        #' initialize(CS50Gram,3) initializes up to trigram
        stopifnot (is.numeric(n), is.finite(n), n > 0, n <= maxlen,
                   is.numeric(level), is.finite(level), level >= 0, is.character(x))
        extract <- magrittr::extract
        filelist<-list()
        filelist <- sapply(1:n, function(i) {filelist[i] <- paste0(x,"Gram",i,"X")})
        x <- list()
        print ("Reading Gram:")
        x <- sapply (seq_along (filelist), function(i) { 
                paste( z, filelist[i], sep="/") %>% paste0(".RData") -> filename
                load( file = filename)
                print (c(filelist[i],", "))
                x[[i]] <- Vocabulary  })  # only retain the Gram
        if (level > 0) { x <- pblapply(1:n, function(i) { x[[i]] %<>% extract(.> level) } )
                         x <- lapply(1:n, function(i) { x[[i]] %<>% extract %>% -level } ) }
        x
}, options = list( optimize = 3) )

log10_ceiling <- function(x) { 10^(ceiling(log10(x))) }
# plyr::round_any(x = 345, accuracy = 1000, f = ceiling) 

Good_Turing <- cmpfun(function(x){ 
        #' Good-Turing applied to a Gram
        #' @param x the Gram to evaluate as a list of named integers
        #' @return a list of bins(r) and their corresponding counts (Zr)
        #' @examples 
        mutate <- dplyr::mutate
        stopifnot (length(x)>0)
        y <- x    # work on a backup
        i <- as.integer(tail( y, 1 ) )
        k <- 1L
        df<-data.frame(r = 0L, Nr = 0L, Zr = 0L, TE = NA, GTE = NA, sdTE = 0 )
        while (i < as.numeric( head( y, 1 ) ) ) { 
                df[k,1] <- i # number of times seen (r)
                df[k,2] <-length(y[y == i]) # Cardinal or Number of items (Nr) in the 'times seen' (r) bin
                y <- y[y > i] 
                k <- k + 1L
                i <- as.integer( tail( y, 1) )
        }
        df[k,1] <- i
        df[k,2] <- length( y[y == i])
        # we now have the bin's data frame
        c0 <- df[1,2] ; p0 <- c0/sum(df[,2])    # estimate of unseen species count and probability
        lstar <- 2*df[2,2]/df[1,2]              # Turing estimator: if >1 the class is closed
        # calculate Zr using Church &Gale proposed method: Zr=Nr/(0.5*(t-q))
        len <- nrow(df)
        df[,3] <- c (df[1,3]<-2*df[1,2]/df[2,1],
                     sapply ( 2L:(len-1L), function(i) {df[i,3]<-2*df[i,2]/(df[(i+1L),1]-df[(i-1L),1])}),
                     df[len,3]<-df[len,2]/(df[len,1]-df[(len-1L),1])) 
        plot(x = df[,1], y = df[,2], type = 'l', xlim = (c(1, 1E7 ) ), ylim = (c(1E-4, 1E6 ) ), log = "xy" ,
                       main = "Good-Turing Method", xlab = "r,Frequency", ylab = "Nr Frequency of Frequency", col = "red" )
        lines(x = df[,1], y = df[,3], type = 'l', xlim = (c(1, 1E7 ) ), ylim = (c(1E-4, 1E6 ) ), col = "green" )
        # now compute the Turing Estimates (only possible for continuous bins r!)
        i <- 1L
        while (df[i,1] == as.integer(rownames(df[i,]))) {df[i,4] <-df[(i+1L),1]*df[(i+1L),2]/df[i,2]; i <- i + 1L}
        # now compute the Good Turing Estimates from log10(Zr)=a+b*Log10(r) fit
        df %<>% mutate(lr = log10(r), lZr = log10(Zr) )
        m <- lm(formula = lZr ~ lr, data = df)
        # compute rstar = r*(1+1/r)^(b+1)
        df$GTE<-df[,1]*(1 + 1/df[,1])^(1 + m$coefficients[2])
        # now compute SD(TE)
        df$sdTE<-sapply(1:len, function(i){ ((1 + df$r[i])/df$Nr[i])*sqrt(df$Nr[(i+1L)]*(1 + df$Nr[(i+1L)]/df$Nr[i]))})
        # now locate the swiching position based on p < 0.05 or p < 0.10 criterion
        criterion <- c(1.96, 1.65) ; swpoint <- 0L ; i <- 1L
        while ((swpoint < 1L) & !(i > length(criterion[i]))) {
                df$sw <- ifelse (abs ( (df$TE - df$GTE )/df$sdTE ) > criterion[i] , TRUE, FALSE )
                swpoint <- max(1L,rownames( head( df[which(df$sw == TRUE),],1)))
                i <- i + 1L
        }
        # we now have a switchpoint to merge TE and GTE
        df$sdTE <- df$lr <- df$lZr <- df$sw <-NULL  # not needed any longer
        df$rstar <- df$GTE
        if (swpoint > 1L) {df$rstar[1L:swpoint-1L]<-df$TE[1L:swpoint-1L]}
        # now calculate the re-normalization factor (1-N1/N) only if swpoint > 1L
        renorm <- ifelse (swpoint == 1L, 1L, 1L-df$Nr[1]/sum(df$Nr))
        # now calculate the Simple Good Turing (SGT) estimate
        df$SGT <- df$rstar*renorm
        ymax <- 10^( ceiling( log10( max( df[,7]/df[,1] ) ) ) )
        plot(x = df[,1], y = df[,7]/df[,1], type = 'l',xlim = (c(1, 1E8) ), ylim = (c(0, 1) ),log = "x",
             main = "Simple Good-Turing Estimate", xlab = "r,Frequency", ylab = "r*/r Relative AdjustedFrequency", col = "red" )
        # now we recompute a integer SGT to use for our grams
        df$RSGT <- as.integer(0.5 + df$SGT)
        df$trial <- df$SGT*df$Zr
        ymax <- 10^( ceiling( log10(max(df[,9]) ) ) )
        plot(x = df[,1], y = df[,9], type = 'l', xlim = (c(1, 1E8) ), ylim = (c(1, ymax) ), log = "xy",
             main = "Simple Good-Turing Estimate N.c*", xlab = "r,Frequency", ylab = "N c* 'constant'", col = "red" )
        df$Zr <- df$TE <- df$GTE <- df$rstar <- df$trial <- NULL
        # now we need to fold back smoothed counts in original Gram
        xs <- rownames( head( df[df$r != df$RSGT,],1))    # easy find of index to start decrementing count
        x[x >= xs] %<>% -1L
        x
}, options = list( optimize = 3) ) 

SGT_Normalize_Gram <- cmpfun(function(x) {
        #' Applies Simple Good Turing algorithm to all grams from Gram list
        #' @param x the gram list
        #' @return a list of SGT corrected grams in x 
        #' @examples 
        #' Normalize_Gram(Gram)   iApplies SGT to all grams in the Gram list
        n <- length(x)
        stopifnot (n > 0)
        for (i in 1:n) { x[[i]] %<>% Good_Turing }
        x
}, options = list( optimize = 3) )

Trim_Gram <- cmpfun(function(x,level) {
        #' Applies Trim level threshold specified to all grams from Gram list
        #' @param x the gram list
        #' @param level the threshold to trrim
        #' @return a list of SGT corrected grams in x 
        #' @examples 
        #' Trim_Gram(Gram)   Applies trim to all grams in the Gram list
        extract <- magrittr::extract
        n <- length(x)
        stopifnot (n > 0L,level > 0L)
        x <- pblapply(1:n, function(i) { x[[i]] %<>% extract(.>level)  } )
        x <- lapply(1:n, function(i) { x[[i]] %<>% extract %>% -level } )
        x
}, options = list( optimize = 3) )

finalize <- cmpfun(function(x,y,z,n=maxlen) {
        #' save('GTCS50GramX0.RData')
        #' @param x the gram list
        #' @param y the name of the grambase (e.g. "CS50") 
        #' @param z the size of the grambase (eg) gramdir, xgramdir, lgramdir, mgramdir, sgramdir, tgramdir, ngramdir
        #' @param n integer the number of grams to save
        #' @return a list of grams in x 
        #' @examples 
        #' finalize(Gram,GTCS50)   finalizes up to 7-gram 
        #' finalize(Gram,GTCS50,3) finalizes up to trigram
        stopifnot (is.integer(n), is.finite(n), n > 0L, n <= maxlen, is.character(y))
        extract <- magrittr::extract
        filelist <- list()
        filelist <- sapply(1:n, function(i) {filelist[i] <- paste0(y, "Gram", i, "0")})
        print ("Writing Gram:")
        for (i in 1:length (filelist)) { 
                paste( z, filelist[i], sep="/") %>% paste0(".RData") -> filename
                x[[i]] -> Vocabulary
                print (c(filelist[i],", "))
                save(list=c("Vocabulary"), file = filename)
        }  # only retain the Gram
        x
}, options = list( optimize = 3) )

finalizeX <- cmpfun(function(x,y,z,n=maxlen) {
        #' save('GTCS50GramX0.RData')
        #' @param x the gram list
        #' @param y the name of the grambase (e.g. "CS50") 
        #' @param z the size of the grambase (eg) gramdir, xgramdir, lgramdir, mgramdir, sgramdir, tgramdir, ngramdir
        #' @param n integer the number of grams to save
        #' @return a list of grams in x 
        #' @examples 
        #' finalize(Gram,GTCS50)   finalizes up to 7-gram 
        #' finalize(Gram,GTCS50,3) finalizes up to trigram
        stopifnot (is.integer(n), is.finite(n), n > 0L, n <= maxlen, is.character(y))
        extract <- magrittr::extract
        filelist <- list()
        filelist <- sapply(1:n, function(i) {filelist[i] <- paste0(y, "Gram", i, "X")})
        print ("Writing Gram:")
        for (i in 1:length (filelist)) { 
                paste( z, filelist[i], sep="/") %>% paste0(".RData") -> filename
                x[[i]] -> Vocabulary
                print (c(filelist[i],", "))
                save(list=c("Vocabulary"), file = filename)
        }  # only retain the Gram
        x
}, options = list( optimize = 3 ))

# Revision 8-21-15
##################################################
#
# Step 1 Simple Good Turing on skip grams
#
##################################################

n <- 6L ; level = 0L  # do not trim any further
location <- data.frame(size = c( "X", "L", "M", "S", "T", "N" ),dirs = c( xgramdir, lgramdir, mgramdir, sgramdir, tgramdir, ngramdir ) )
sapply(1:nrow(location), function(j){
        Gram <- list()
        Gram <- initializeX(paste0(location[j,1],"CS50"),location[j,2], n = maxlen, level = 0L)     # work with original Sampled 50% data up to heptigram
        Gram %<>% SGT_Normalize_Gram
        Gram <- sapply (1:maxlen, function(i) { Gram[[i]] <- setNames( as.integer( Gram[[i]] ), names( Gram[[i]] ) ) } )
        Gram %<>% finalizeX(paste0(location[j,1],"GTCS50"), location[j,2], n = maxlen)                 # save SGT and leveled Samples 80% data up to heptigram
        saveRDS(Gram,file = paste0(location[j,1],"GTCS5X.RDS"))
        })
##################################################
#
# Step 2 Simple Good Turing on full grams
#
##################################################
n <- 6L ; level = 0L  # do not trim any further
location <- data.frame(size = c( "X", "L", "M", "S", "T", "N" ),dirs = c( xgramdir, lgramdir, mgramdir, sgramdir, tgramdir, ngramdir ) )
sapply(1:nrow(location), function(j){
        Gram <- list()
        Gram <- initialize(paste0(location[j,1],"CS50"),location[j,2], n = maxlen, level = 0L)     # work with original Sampled 50% data up to heptigram
        Gram %<>% SGT_Normalize_Gram
        Gram <- sapply (1:maxlen, function(i) { Gram[[i]] <- setNames( as.integer( Gram[[i]] ), names( Gram[[i]] ) ) } )
        Gram %<>% finalize(paste0(location[j,1],"GTCS50"), location[j,2], n = maxlen)                 # save SGT and leveled Samples 80% data up to heptigram
        saveRDS(Gram,file = paste0(location[j,1],"GTCS50.RDS"))
        })

