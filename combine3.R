# combine3.R combines grams and skip-grams
options(mc.cores=4)
memory.limit(size=28000)
setwd("V:/") 
Sys.info()[1:5]                         # system info but exclude login and user info
userdir<-getwd()
# user-defined startup Capstone directory
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

initialize <- cmpfun (function(x,k,n,level=0L) {
        #' load('CS50Gramkn.RData')
        #' @param x the name of the grambase (e.g. "CS80")
        #' @param k the type of gram base (2 for bigrams, 3 for trigrams...) 
        #' @param n integer the number of skip grams to load
        #' @param level integer to trim
        #' @return a list of grams in x 
        #' @examples 
        #' initialize(CS50Gram)   initializes up to 7-gram 
        #' initialize(CS50Gram,3) initializes up to trigram
        stopifnot (is.integer(n), is.finite(n), n > -1L, n <= maxlen,
                   is.integer(level), is.finite(level), level > -1L, 
                   is.integer(k), is.finite(k), k > 1L, is.character(x))
        extract <- magrittr::extract
        n <- n + 1L  # we have the series 0..n to handle
        filelist<-list()
        filelist <- sapply(1:n, function(i) {filelist[i] <- paste0(x,"Gram",k,(i-1L)) } )
        x <- list()
        print ("Reading Gram:")
        x <- sapply (seq_along (filelist), function(i) { 
                paste( gramdir, filelist[i], sep="/") %>% paste0(".RData") -> filename
                print (c(filelist[i],", "))
                load(file=filename)
                x[[i]] <- Vocabulary })
        x <- pblapply(1:n, function(i) { x[[i]] %<>% extract(.> level) })
        x <- lapply(1:n, function(i) { x[[i]] %<>% extract %>% -level })
        x
}, options = list(optimize=3))

merge_grams <- cmpfun (function(x) {
        #' merge ('CS50Gramkn.RData')
        #' @param x the name of the grambase (e.g. "CS80")
        #' @return a list of grams in x 
        #' @examples 
        #' merge_grams(Gram)   merges all grams loaded in the list
        extract <- magrittr::extract
        stopifnot (length(x) > 0)
        n <- length(x)
        y <- x[[1]]
        for (i in 1:n) { y <- c( y, x[[i]] ) } 
        # cleanup repetitions of words of 2 characters or more ... 
        names(y) %<>% gsub("(?:\\b(\\w+)\\b) (?:\\1(?: |$))+",'*', ., ignore.case = TRUE, perl = TRUE, useBytes = TRUE)
        if (n > 2) {names(y) %<>% gsub("^(\\*.*)|(.*\\*.*)|(.*\\*)$",'*', ., ignore.case = TRUE, perl = TRUE, useBytes = TRUE)}
        y %>% aggregate( by = list( names(.) ), FUN = sum, simplify = TRUE ) -> x
        x$x <- as.integer( x$x )
        y <- setNames(as.integer( x$x ), as.list( x$Group.1 ) )
        y %>% sort( decreasing = TRUE ) -> x
        x[names(x) == '*'] <- 0L        # clear the value for the special pattern *
        x <- x[x > 0L]
        x
}, options = list( optimize = 3) )

#################################################
#
# step 0. Initialize
#
#################################################
Gram <- list() 
#################################################
#
# step 1. Unigram
#
#################################################
filename<-paste(gramdir, "CS50Gram10.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 levels for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram1X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram1X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram1X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram1X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram1X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram1X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 2. Bi skipgrams
#
#################################################
Gram <- initialize( "CS50", k = 2L, n = 4L, level = 4L)        # work with Raw Sampled 50% bigrams and skipgrams, trim by 4 levels
Gram %>% merge_grams -> Vocabulary
# save combined grams
filename <-paste(xgramdir,"XCS50Gram2X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram2X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram2X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram2X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram2X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram2X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 3. Tri skipgrams
#
#################################################
Gram <- initialize( "CS50", k = 3L, n = 3L, level = 4L)        # work with Raw Sampled 50% trigrams and skipgrams
Gram %>% merge_grams -> Vocabulary
# save combined grams
filename <-paste(xgramdir,"XCS50Gram3X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram3X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram3X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram3X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram3X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram3X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 4. Quadri skipgrams
#
#################################################
Gram <- initialize( "CS50", k = 4L, n = 2L, level = 4L)        # work with Raw Sampled 50% quadrigrams and skipgrams
Gram %>% merge_grams -> Vocabulary
# save combined grams
filename <-paste(xgramdir,"XCS50Gram4X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram4X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram4X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram4X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram4X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram4X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 5. Quinti skipgrams
#
#################################################
Gram <- initialize( "CS50", k = 5L, n = 1L, level = 4L)        # work with Raw Sampled 50% quintigrams and skipgrams
Gram %>% merge_grams -> Vocabulary
# save combined grams
filename <-paste(xgramdir,"XCS50Gram5X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram5X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram5X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram5X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram5X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram5X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 6. Sexti skipgrams
#
#################################################
filename<-paste(gramdir, "CS50Gram60.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 more level for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram6X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gramVocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 more level for a total of 4 level
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram6X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram6X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram6X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram6X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram6X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 7. Septigram
#
#################################################
#filename<-paste(gramdir, "CS50Gram70.RData", sep = "/" ) ; load( file = filename )
#Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 more level for a total of 4 level
#Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram7X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
#Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
#Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram7X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
#Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
#Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram7X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
#Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
#Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram7X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
#Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
#Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram7X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
#Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
#Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram7X.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# Now handle the regular grams
#
#################################################
#
# step 1. Unigram
#
#################################################
filename<-paste(gramdir, "CS50Gram10.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 levels for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram10.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram10.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram10.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram10.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram10.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram10.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 2. Bigrams
#
#################################################
filename<-paste(gramdir, "CS50Gram20.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 levels for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram20.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram20.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram20.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram20.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram20.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram20.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 3. Tri skipgrams
#
#################################################
filename<-paste(gramdir, "CS50Gram30.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 levels for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram30.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram30.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram30.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram30.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram30.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram30.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 4. Quadri skipgrams
#
#################################################
filename<-paste(gramdir, "CS50Gram40.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 levels for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram40.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram40.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram40.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram40.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram40.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram40.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 5. Quinti skipgrams
#
#################################################
filename<-paste(gramdir, "CS50Gram50.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 levels for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram50.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram50.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram50.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram50.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram50.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram50.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 6. Sexti skipgrams
#
#################################################
filename<-paste(gramdir, "CS50Gram60.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 levels for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram60.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram60.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram60.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram60.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram60.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram60.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
#
# step 7. Septigram
#
#################################################
filename<-paste(gramdir, "CS50Gram70.RData", sep = "/" ) ; load( file = filename )
Vocabulary <- Vocabulary[ Vocabulary > 4L]                      # trim by 4 more level for a total of 4 level
Vocabulary <- Vocabulary - 4L ; filename <-paste(xgramdir,"XCS50Gram70.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Xlarge gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 6 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(lgramdir,"LCS50Gram70.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Large gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 8 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(mgramdir,"MCS50Gram70.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Medium gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 10 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(sgramdir,"SCS50Gram70.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Small gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 2 more levels for a total of 12 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(tgramdir,"TCS50Gram70.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Tiny gram
Vocabulary <- Vocabulary[ Vocabulary > 2L]                      # trim by 1 more levels for a total of 14 levels
Vocabulary <- Vocabulary - 2L ; filename <-paste(ngramdir,"NCS50Gram70.RData",sep="/") ; save(list =c ( "Vocabulary" ), file = filename) # save as a combined Nano gram
#################################################
