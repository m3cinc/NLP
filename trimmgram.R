# trimmgram.R
options(mc.cores=4)
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

initialize <- function(x,n) {
        #' load('CS80GramX0.RData')
        #' @param x the name of the grambase (e.g. "CS80") 
        #' @param n integer the number of grams to load
        #' @return a list of grams in x 
        #' @examples 
        #' initialize(CS80Gram)   initializes up to 7-gram 
        #' initialize(CS80Gram,3) initializes up to trigram
        stopifnot (is.numeric(n), is.finite(n), n > 0, n <= maxlen, is.character(x))
        extract <- magrittr::extract
        filelist<-list()
        filelist <- sapply(1:n, function(i) {filelist[i] <- paste0(x,"Gram",i,"0")})
        x <- list()
        x <- sapply (seq_along (filelist), function(i) { 
                paste( gramdir, filelist[i], sep="/") %>% paste0(".RData") -> filename
                load(file=filename)
                x[[i]] <- Vocabulary  })  # only retain the Gram
#        x <- sapply(1:n, function(i) { x[[i]] %<>% unlist %>% extract(.>5) } )
#        x <- sapply(1:n, function(i) { x[[i]] %<>% unlist %>% extract %>% -5L } )
        x
}
initialize_Level3<-cmpfun(initialize,options = list(optimize=3))

Good_Turing <-function(x){ 
        #' Good-Turing applied to a Gram
        #' @param x the Gram to evaluate as a named integer list
        #' @return a list of bins(r) and their corresponding counts (Zr)
        #' @examples 
        #' Good_Turing(x)   perform algorithm on list x
        mutate <- dplyr::mutate
        stopifnot (length(x)>0)
        y <- x    # work on a backup
        i <- as.numeric(tail(y,1))
        k <- 1
        df<-data.frame(r=0L,Nr=0L,Zr=0,TE=NA,GTE=NA,sdTE=0)
        while (i < as.numeric(head(y,1))) { 
                df[k,1] <- i # number of times seen (r)
                df[k,2] <-length(y[y==i]) # Cardinal or Number of items (Nr) in the 'times seen' (r) bin
                y <- y[y > i] 
                k<-k+1
                i <- as.numeric(tail(y,1))
        }
        df[k,1]<-i
        df[k,2]<-length(y[y==i])
        # we now have the bin's data frame
        c0 <- df[1,2] ; p0 <- c0/sum(df[,2])    # estimate of unseen species count and probability
        lstar <- 2*df[2,2]/df[1,2]              # Turing estimator: if >1 the class is closed
        # calculate Zr using Church &Gale proposed method: Zr=Nr/(0.5*(t-q))
        len <- nrow(df)
        df[,3] <- c (df[1,3]<-2*df[1,2]/df[2,1],
                     sapply ( 2L:(len-1L), function(i) {df[i,3]<-2*df[i,2]/(df[i+1,1]-df[i-1,1])}),
                     df[len,3]<-df[len,2]/(df[len,1]-df[(len-1L),1])) 
        plot(x=df[,1],y=df[,2], type='l', xlim=(c(1,1E7)), ylim=(c(1E-4,1E6)), log = "xy" ,
                       main="Good-Turing Method",xlab="r,Frequency",ylab="Nr Frequency of Frequency",col="red")
        lines(x=df[,1],y=df[,3], type='l',xlim=(c(1,1E7)),ylim=(c(1E-4,1E6)),col="green")
        # now compute the Turing Estimates (only possible for continuous bins r!)
        i<-1
        while (df[i,1] == as.numeric(rownames(df[i,]))) {df[i,4] <-df[i+1L,1]*df[i+1L,2]/df[i,2]; i<-i+1}
        # now compute the Good Turing Estimates from log10(Zr)=a+b*Log10(r) fit
        df %<>% mutate(lr=log10(r),lZr=log10(Zr))
        m<-lm(formula = lZr ~ lr, data=df,  method="qr")
        summary(m)
        # compute rstar = r*(1+1/r)^(b+1)
        df$GTE<-df[,1]*(1+1/df[,1])^(1+m$coefficients[2])
        # now compute SD(TE)
        df$sdTE<-sapply(1:len, function(i){ ((1+df$r[i])/df$Nr[i])*sqrt(df$Nr[i+1]*(1+df$Nr[i+1]/df$Nr[i]))})
        # now locate the swiching position based on p < 0.05 criterion
        criterion <- c(1.96,1.65) ; swpoint<-0L; i<-1L
        while ((swpoint < 1L) & !(i > length(criterion[i]))) {
                df$sw <- ifelse (abs((df$TE-df$GTE)/df$sdTE) > criterion[i] , TRUE, FALSE )
                swpoint <-rownames(head(df[which(df$sw == TRUE),],1))
                i<-i+1
        }
        # we now have a switchpoint to merge TE and GTE
        df$sdTE <- df$lr <- df$lZr <- df$sw <-NULL  # not needed any longer
        df$rstar<-df$GTE
        df$rstar[1:swpoint-1L]<-df$TE[1:swpoint-1L]
        # now calculate the re-normalization factor (1-N1/N)
        renorm <- 1L-df$Nr[1]/sum(df$Nr)
        # now calculate the Simple Good Turing (SGT) estimate
        df$SGT<-df$rstar*renorm
        plot(x=df[,1],y=df[,7]/df[,1], type='l',xlim=(c(1,1E8)),ylim=(c(0.5,1)),log = "x",
             main="Simple Good-Turing Estimate",xlab="r,Frequency",ylab="r*/r Relative AdjustedFrequency",col="red")
        # now we recompute a integer SGT to use for our grams
        df$RSGT<-ceiling(df$SGT)
        df$trial<-df$SGT*df$Zr
        plot(x=df[,1],y=df[,9], type='l',xlim=(c(1,1E8)),ylim=(c(1,50000)),log = "xy",
             main="Simple Good-Turing Estimate N.c*",xlab="r,Frequency",ylab="N c* 'constant'",col="red")
        df$Zr <- df$TE <- df$GTE <- df$rstar <- df$trial <- NULL
        df$RSGT <- as.integer(df$RSGT)
        # now we need to fold back smoothed counts in original Gram
        x <- pbsapply(len:1, function(i){ setNames(rep.int(df[i,4],length(x[x==df[i,1]])), names(x[x==df[i,1]])) } )
}
Good_Turing_Level3<-cmpfun(Good_Turing,options = list(optimize=3))

SGT_Normalize_Gram <- function(x) {
        #' Applies Simple Good Turing algorithm to all grams from Gram list
        #' @param x the name of the gram list
        #' @return a list of SGT corrected grams in x 
        #' @examples 
        #' Normalize_Gram(Gram)   iApplies SGT to all grams in the Gram list
        n <- length(x)
        stopifnot (n > 0)
        x <- pbsapply(1:n, function(i) { print(c("Applying SGT at Gram ",i,"..."))
                                         x[i] %<>% unlist %>% Good_Turing_Level3 } )
}
SGT_Normalize_Gram_Level3<-cmpfun(SGT_Normalize_Gram,options = list(optimize=3))

Trim_Gram <- function(x,level) {
        #' Applies Trim level threshold specified to all grams from Gram list
        #' @param x the name of the gram list
        #' @param level the threshold to trrim
        #' @return a list of SGT corrected grams in x 
        #' @examples 
        #' Trim_Gram(Gram)   Applies trim to all grams in the Gram list
        extract <- magrittr::extract
        n <- length(x)
        stopifnot (n > 0,level > 0)
        x <- pbsapply(1:n, function(i) { x[i] %<>% unlist %>% extract(.>level)  } )
}
Trim_Gram_Level3<-cmpfun(Trim_Gram,options = list(optimize=3))

finalize <- function(x,n) {
        #' save('GTCS80GramX0.RData')
        #' @param x the name of the grambase (e.g. "CS80") 
        #' @param n integer the number of grams to load
        #' @return a list of grams in x 
        #' @examples 
        #' finalize(GTCS80Gram)   finalizes up to 7-gram 
        #' finalize(GTCS80Gram,3) finalizes up to trigram
        stopifnot (is.numeric(n), is.finite(n), n > 0, n <= maxlen, is.character(x))
        extract <- magrittr::extract
        filelist<-list()
        filelist <- sapply(1:n, function(i) {filelist[i] <- paste0(x,"Gram",i,"0")})
        for (i in 1:length (filelist)) { 
                paste( gramdir, filelist[i], sep="/") %>% paste0(".RData") -> filename
                x[[i]] -> Vocabulary
                save(list=c("Vocabulary"),file=filename)
        }  # only retain the Gram
        x
}
finalize_Level3<-cmpfun(initialize,options = list(optimize=3))


n<-4L;level=5L
Gram <- list()
Gram <- initialize_Level3("S80",n=4L)       # work with Raw Sampled 80% data up to quadrigram
Gram <- SGT_Normalize_Gram_Level3(Gram)
Gram <- Trim_Gram_Level3(Gram,level=5L)
Gram <- finalize_Level3("GTCS80",n=4L)      # save SGT and leveled Samples 80% data up to quadrigram
