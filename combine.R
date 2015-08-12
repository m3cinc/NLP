# combine.R combines grams and skip-gramsoptions(mc.cores=4)
memory.limit(size=28000)
setwd("W:/") 
Sys.info()[1:5]                         # system info but exclude login and user info
userdir<-getwd()
# user-defined startup Capstone directory
library(dplyr)
library(magrittr)
library(stringi)
library(pbapply)
library(rbenchmark)
library(microbenchmark)
library(proftools)
library(compiler)
local <- TRUE                           # to avoid data redownload
figdir <- "./DSC_Milestone_files/figure-html"
datadir <- "./data"
datasetsdir <- paste0(datadir,"/datasets")
userdatadir <- paste0(datadir,"/Final/en_US")
cleandir <- paste0(datadir,"/clean")
sampledir <- paste0(datadir,"/sample")
process_fulldir <- paste0(datadir,"/processed")
textdir <- paste0(datadir,"/text")
gramdir <- paste0(datadir,"/gram")
maxlen <-7
if (!file.exists(datadir)) { dir.create(datadir) }  #  data will reside in datadir
if (!file.exists(datasetsdir)) { dir.create(datasetsdir) }  #  datasets will reside in datasetsdir
if (!file.exists(cleandir)) { dir.create(cleandir) }  # clean will reside in cleandir
if (!file.exists(sampledir)) { dir.create(sampledir) }  # sample will reside in sampledir
if (!file.exists(process_fulldir)) { dir.create(process_fulldir) }  # process_full will reside in process_fulldir
if (!file.exists(textdir)) { dir.create(textdir) }  #  data will reside in datadir
if (!file.exists(gramdir)) { dir.create(gramdir) }  #  gram will reside in datadir

initialize <- function(x,k,n,level=0) {
        #' load('CS80Gramkn.RData')
        #' @param x the name of the grambase (e.g. "CS80")
        #' @param m the type of gram base (2 for bigrams, 3 for trigrams...) 
        #' @param n integer the number of skip grams to load
        #' @param level integer to trim
        #' @return a list of grams in x 
        #' @examples 
        #' initialize(CS80Gram)   initializes up to 7-gram 
        #' initialize(CS80Gram,3) initializes up to trigram
        stopifnot (is.numeric(n), is.finite(n), n >= 0, n <= maxlen,
                   is.numeric(level), is.finite(level), level >= 0, 
                   is.numeric(k), is.finite(k), k > 1, is.character(x))
        extract <- magrittr::extract
        n<-n+1  # we have the series 0..n to handle
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
} 
initialize_Level3<-cmpfun(initialize,options = list(optimize=3))

log10_ceiling <- function(x) { 10^(ceiling(log10(x))) }
# plyr::round_any(x = 345, accuracy = 1000, f = ceiling) 

merge_grams <-function(x) {
        #' merge ('CS80Gramkn.RData')
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
        names(y) %<>% gsub("(?:\\b(\\w+)\\b) (?:\\1(?: |$))+",'*',.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE)
        if (n > 2) {names(y) %<>% gsub("^(\\*.*)|(.*\\*.*)|(.*\\*)$",'*',.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE)}
        y %>% aggregate(by=list(names(.)),FUN=sum,simplify=TRUE) -> x
        x$x<-as.integer(x$x)
        y <- setNames(as.integer(x$x),as.list(x$Group.1))
        # clear the value for the special pattern *
        y %>% sort(decreasing=TRUE) -> x
        x[names(x)=='*'] <- 0L
        X <- x[x>0]
        x
} 
merge_grams_Level3<-cmpfun(merge_grams,options = list(optimize=3)) 

n<-4L ; level=8L ; k=2
Gram <- list() 

filename<-paste(gramdir,"CS80Gram10.RData",sep="/") 
load(file=filename)
x<-Vocabulary
x<-x[x>8L]
x<-x-8L
x %>% aggregate(by=list(names(.)),FUN=sum,simplify=TRUE) -> x
x$x<-as.integer(x$x)
x <- setNames(x$x,as.list(x$Group.1))
# clear the value for the special pattern *
x %<>% sort(decreasing=TRUE) 
x[names(x)=='*']<-0
x[x>0]->x
Vocabulary<-x
# save combined grams
filename <-paste(gramdir,"CS80Gram1X.RData",sep="/") 
save(list=c("Vocabulary"),file=filename)

Gram <- initialize_Level3("CS80",k=2,n=4L,level=8L)        # work with Raw Sampled 80% bigrams and skipgrams
Gram %>% merge_grams_Level3 -> Vocabulary
# save combined grams
filename <-paste(gramdir,"CS80Gram2X.RData",sep="/") 
save(list=c("Vocabulary"),file=filename)

Gram <- initialize_Level3("CS80",k=3,n=1L,level=8L) 
Gram %>% merge_grams_Level3 -> Vocabulary
# save combined grams
filename <-paste(gramdir,"CS80Gram3X.RData",sep="/") 
save(list=c("Vocabulary"),file=filename)

filename<-paste(gramdir,"CS80Gram40.RData",sep="/") 
load(file=filename)
x<-Vocabulary
x<-x[x>8L]
x<-x-8L
# cleanup repetitions of words of 2 characters or more ... 
names(x) %<>% gsub("(?:\\b(\\w+)\\b) (?:\\1(?: |$))+",'*',.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE)
names(x) %<>% gsub("^(\\*.*)|(.*\\*.*)|(.*\\*)$",'*',.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE)
x %>% aggregate(by=list(names(.)),FUN=sum,simplify=TRUE) -> x
x$x<-as.integer(x$x)
x <- setNames(x$x,as.list(x$Group.1))
# clear the value for the special pattern *
x %<>% sort(decreasing=TRUE) 
x[names(x)=='*']<-0
x[x>0]->x
Vocabulary<-x
# save combined grams
filename <-paste(gramdir,"CS80Gram4X.RData",sep="/") 
save(list=c("Vocabulary"),file=filename)



