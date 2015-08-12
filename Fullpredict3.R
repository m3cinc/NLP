# fullpredict3.R
options(mc.cores=4)
memory.limit(size=28000)
setwd("W:/") 
Sys.info()[1:5]                         # system info but exclude login and user info
userdir<-getwd()                        # user-defined startup Capstone directory
library(magrittr)
library(stringi)
library(tm)
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

library(rbenchmark)
library(microbenchmark)
library(proftools)
library(magrittr)
library(stringi)
library(compiler)


initialize <- function(x,n,level=0L) {
        #' load('CS80GramX0.RData')
        #' @param x the name of the grambase (e.g. "CS80") 
        #' @param n integer the number of grams to load
        #' @param level integer to trim
        #' @return a list of grams in x 
        #' @examples 
        #' initialize(CS80Gram)   initializes up to 7-gram 
        #' initialize(CS80Gram,3) initializes up to trigram
        stopifnot (is.numeric(n), is.finite(n), n > 0, n <= maxlen,
                   is.numeric(level), is.finite(level), level >= 0, is.character(x))
        extract <- magrittr::extract
        filelist<-list()
        filelist <- sapply(1:n, function(i) {filelist[i] <- paste0(x,"Gram",i,"X")})
        x <- list()
        print ("Reading Gram:")
        x <- sapply (seq_along (filelist), function(i) { 
                paste( gramdir, filelist[i], sep="/") %>% paste0(".RData") -> filename
                load(file=filename)
                print (c(filelist[i],", "))
                x[[i]] <- Vocabulary  })  # only retain the Gram
        x <- pblapply(1:n, function(i) { x[[i]] %<>% extract(.> level) } )
        x <- lapply(1:n, function(i) { x[[i]] %<>% extract %>% -level } )
        x
}
initialize_Level3<-cmpfun(initialize,options = list(optimize=3))

prepare <- function(x)  {
        #' initialize
        #' @param x the input string
        #' @return a pattern list of regex encoded tokens in x 
        #' @examples 
        #' initialize(mystring)   initializes, tokenizes and Regex prepares a patern list for grep 
        ## Prepare x and return pattern list in x
        stopifnot(is.character(x), is.finite(n), n > 0, n <= maxlen, is.character(x))
        #' To avoid :: calls 
        stri_split_boundaries <- stringi::stri_split_boundaries 
        stri_join <- stringi::stri_join  
        stri_flatten <- stringi::stri_flatten
        options <- stringi::stri_opts_brkiter( type="word", skip_word_none = TRUE, skip_word_number = FALSE )
        # clean input
        x %<>%  gsub("[^A-Za-z\' ]","",.) %>%  tolower %>%
                stri_split_boundaries(., opts_brkiter=options) %>% unlist
        x %>% length -> len
        # Split into word tokens
        x <- sapply (1L:len, function(i) {stri_flatten(tail(x,i),collapse=" ")})
        x %<>% gsub("\\'","\\\\\\'",.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE) 
        x <- sapply (1:(len+1L), function(i) {
                x[i] <- ifelse (i > 1L, paste0("^(\\b",x[i-1L],")(?=[[:space:]])(.*)"),paste0("^(\\b",x[1],")(?![[:punct:]]|[[:alnum:]])") ) } )
} 
prepare_Level3<-cmpfun(prepare,options = list(optimize=3))


suggest_nlpmodel <- function(x,model=1) {
        #' natural language processor to suggest continuation word
        #' @param x the input string in  regex encoded tokens in x
        #' @param model defines the model type as follows
        #'        model = 1 for 'Simple Backoff' with s = 0.4 using full-grams
        #'        model = 2 for 'weighted'Compounded' with backoff using full grams
        #'        model = 3 for 'Simple Backoff' with s = 0.4 using skip-grams
        #' @return a list of suggestion strings x in decreasing scores
        #' @examples 
        #' initialize(mystring)   initializes, tokenizes and Regex prepares a patern list for grep 
        ## Prepare x and return pattern list in x
        pattern <- x
        n <- length(x)
        x <- list()
        y <- list()
        x <- sapply (1:(n), function(i) {
                .subset2(Gram,i)[grep(pattern[i],names(Gram[[i]]),ignore.case=TRUE,perl=TRUE,value=FALSE,useBytes=TRUE)] } )
        # use default 'Stupid Backoff' factor
        s <- 0.4
        ## Predict Algorithm
        # initialize
        x <- Filter(length,x)
        n <- length(x)
        x <- sapply(1:(n), function(i) { if (i > 1) {x[[i]] <- s^(n-i)*x[[i]]/sum(x[[i]])} else {x[[i]] <- s^(n-i)*x[[i]]/sum(Gram[[i]])} } ) 
        x <- sapply(1:(n), function(i) { x[[i]] %<>% head(10)})         # only work on top of the lists
        y <- sapply(1:(n), function(i) { y[[i]] <- sub(pattern[i],"\\2",names(x[[i]]),ignore.case=TRUE,perl=TRUE,useBytes=TRUE) } ) 
        y <- sapply(1:(n), function(i) { gsub("(^ )(.*)","\\2",y[[i]],ignore.case=TRUE,perl=TRUE,useBytes=TRUE) } )
        y[[1]] <- "?"    # just in case we couldnot match anything
        for (i in 1:(n)) { names(x[[i]]) <- y[[i]] }
        x %<>% rev %>% unlist 
        switch ( as.character(model) ,
        "1" = { },
        "2" = { x %>% data.frame(suggestion=names(.),score=.,stringsAsFactors=FALSE) -> z
                z <- aggregate(z$score,list(z$suggestion),sum)
                names (z) <- c("suggestion","score")
                x <- z[,2] ; names(x) <- z[,1] ; rm(z,n)
                x %<>% sort(decreasing=TRUE) } )
        x %<>% head(5) %>% names  
}         
suggest_nlpmodel_Level3<-cmpfun(suggest_nlpmodel,options = list(optimize=3))

google_predictor <- function(x) {
        #' retrieve suggestion from suggestqueries.google.com
        #' @param x the input string (must be char and of length >0)
        #' @return a list of suggestion strings in x in decreasing score rank 
        #' @examples 
        #' google_predictor('you hadn't time to take a') retruns list x: [1] "class"   "bath"    "break"   "nap"     "picture"
        ## Prepare x and return pattern list in x
        stopifnot (is.character(x), nchar(x)>0)        
        x %<>% tolower %>% gsub("[^a-z\\' ]","",.) %>% strsplit(.,' ') %>% unlist %>% tail(3) %>% paste(.,collapse = ' ')
        x %>% strsplit(.,' ') %>% unlist  %>% length -> offset
        x %>% URLencode %>% paste0('http://suggestqueries.google.com/complete/search?client=firefox&q=',.,'%20') -> google_url
        toString(sample(100000:999999, 1)) -> random_destination
        suppressWarnings(download.file(google_url,random_destination,'internal',quiet=TRUE,mode="w",cacheOK=FALSE))
        readChar(random_destination,file.info(random_destination)$size) -> x
        unlink (random_destination) ;rm(random_destination,google_url)
        x %>% gsub('((u0027))',"\\'",.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE) %>%
               gsub('(\\\\)','',.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE) %>%
               gsub('(( *[[:punct:]]{2,}))','\\\r',.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE) %>%
               strsplit(.,split='\r',fixed=TRUE) %>% unlist %>% extract (2:length(.)) -> x
        x %>% extract (1) %>% paste0('(.*)(',.,')( +)(.*)') -> pattern
        if (length(x)==1) {x <- "?" ; return(x)}
        x %>% extract(2:length(.)) %>% gsub(pattern,'\\4',.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE) %>%
        gsub('(^\\w+)(.*)','\\1',.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE) %>% unlist %>% table %>% sort(decreasing=TRUE) %>%
        head (5) %>% names -> x
        return (x)
} 

google_predictor_Level3<-cmpfun(google_predictor,options = list(optimize=3))

######### Example follow
z<-runif(1e3)
Rprof("Rprof.out")
Gram<-initialize_Level3("GTCS80",n=4,level=0)      # work with Cleaned Sampled 80% data up to trigram
z<-runif(1e3)
Rprof(NULL)
summaryRprof("Rprof.out")

Rprof("Rprof.out")
x<-"live and I'd"
n<-sapply(strsplit(x, " "), length)     # count words
x %>% prepare_Level3 %>% suggest_nlpmodel_Level3(model=1) -> m1;m1
x %>% google_predictor -> m2;m2
Rprof(NULL)
summaryRprof("Rprof.out")

Rprof("Rprof.out")
x<-"me about his"
x %>% google_predictor -> m2;m2
Rprof(NULL)
summaryRprof("Rprof.out")

x<-"for you I'd live and I'd"
pattern <-prepare(x)
pattern
# [1] "^(\\bi\\'d)(?![[:punct:]]|[[:alnum:]])"                    
# [2] "^(\\band i\\'d)(?![[:punct:]]|[[:alnum:]])"                
# [3] "^(\\blive and i\\'d)(?![[:punct:]]|[[:alnum:]])"         
# [4] "^(\\bid live and i\\'d)(?![[:punct:]]|[[:alnum:]])"        
# [5] "^(\\byou id live and i\\'d)(?![[:punct:]]|[[:alnum:]])"    
# [6] "^(\\bfor you id live and i\\'d)(?![[:punct:]]|[[:alnum:]])"

x<-" he started telling me about his"
pattern <-prepare(x)
pattern
# [1] "^(\\bstarted)(?![[:punct:]]|[[:alnum:]])"                           
# [2] "^(\\bhe started)(?![[:punct:]]|[[:alnum:]])"                        
# [3] "^(\\band he started)(?![[:punct:]]|[[:alnum:]])"                    
# [4] "^(\\bdessert and he started)(?![[:punct:]]|[[:alnum:]])"            
# [5] "^(\\babout dessert and he started)(?![[:punct:]]|[[:alnum:]])"      
# [6] "^(\\basked about dessert and he started)(?![[:punct:]]|[[:alnum:]])"

x<-"anything to see arctic monkeys this"
pattern <-prepare(x)
pattern
# [1] "^(\\bthis)(?![[:punct:]]|[[:alnum:]])"                               
# [2] "^(\\bmonkeys this)(?![[:punct:]]|[[:alnum:]])"                       
# [3] "^(\\barctic monkeys this)(?![[:punct:]]|[[:alnum:]])"                
# [4] "^(\\bsee arctic monkeys this)(?![[:punct:]]|[[:alnum:]])"            
# [5] "^(\\bto see arctic monkeys this)(?![[:punct:]]|[[:alnum:]])"         
# [6] "^(\\banything to see arctic monkeys this)(?![[:punct:]]|[[:alnum:]])"

x<-"a hug and helps reduce your"
pattern <-prepare(x)
pattern
# [1] "^(\\byour)(?![[:punct:]]|[[:alnum:]])"                       
# [2] "^(\\breduce your)(?![[:punct:]]|[[:alnum:]])"                
# [3] "^(\\bhelps reduce your)(?![[:punct:]]|[[:alnum:]])"          
# [4] "^(\\band helps reduce your)(?![[:punct:]]|[[:alnum:]])"      
# [5] "^(\\bhug and helps reduce your)(?![[:punct:]]|[[:alnum:]])"  
# [6] "^(\\ba hug and helps reduce your)(?![[:punct:]]|[[:alnum:]])"


x<-"you hadn't time to take a"
pattern <-prepare(x)
pattern
# [1] "^(\\ba)(?![[:punct:]]|[[:alnum:]])"                          
# [2] "^(\\btake a)(?![[:punct:]]|[[:alnum:]])"                     
# [3] "^(\\bto take a)(?![[:punct:]]|[[:alnum:]])"                  
# [4] "^(\\btime to take a)(?![[:punct:]]|[[:alnum:]])"             
# [5] "^(\\bhadn\\'t time to take a)(?![[:punct:]]|[[:alnum:]])"    
# [6] "^(\\byou hadn\\'t time to take a)(?![[:punct:]]|[[:alnum:]])"

x<-"and a jury to settle the"
pattern <-prepare(x)
pattern
# [1] "^(\\bthe)(?![[:punct:]]|[[:alnum:]])"                              
# [2] "^(\\bsettle the)(?![[:punct:]]|[[:alnum:]])"                       
# [3] "^(\\bto settle the)(?![[:punct:]]|[[:alnum:]])"                    
# [4] "^(\\bjury to settle the)(?![[:punct:]]|[[:alnum:]])"               
# [5] "^(\\ba jury to settle the)(?![[:punct:]]|[[:alnum:]])"             
# [6] "^(\\band a jury to settle the)(?![[:punct:]]|[[:alnum:]])"         

x<-"of bags of groceries in each"
pattern <-prepare(x)
pattern
# [1] "^(\\beach)(?![[:punct:]]|[[:alnum:]])"                        
# [2] "^(\\bin each)(?![[:punct:]]|[[:alnum:]])"                     
# [3] "^(\\bgroceries in each)(?![[:punct:]]|[[:alnum:]])"           
# [4] "^(\\bof groceries in each)(?![[:punct:]]|[[:alnum:]])"        
# [5] "^(\\bbags of groceries in each)(?![[:punct:]]|[[:alnum:]])"   
# [6] "^(\\bof bags of groceries in each)(?![[:punct:]]|[[:alnum:]])"

x<-"perfect from the bottom to the"
pattern <-prepare(x)
pattern
# [1] "^(\\bthe)(?![[:punct:]]|[[:alnum:]])"                           
# [2] "^(\\bto the)(?![[:punct:]]|[[:alnum:]])"                        
# [3] "^(\\bbottom to the)(?![[:punct:]]|[[:alnum:]])"                 
# [4] "^(\\bthe bottom to the)(?![[:punct:]]|[[:alnum:]])"             
# [5] "^(\\bfrom the bottom to the)(?![[:punct:]]|[[:alnum:]])"        
# [6] "^(\\bperfect from the bottom to the)(?![[:punct:]]|[[:alnum:]])"

x<-"with imagination and bruises from playing"
pattern <-prepare(x)
pattern
# [1] "^(\\bplaying)(?![[:punct:]]|[[:alnum:]])"                                  
# [2] "^(\\bfrom playing)(?![[:punct:]]|[[:alnum:]])"                             
# [3] "^(\\bbruises from playing)(?![[:punct:]]|[[:alnum:]])"                     
# [4] "^(\\band bruises from playing)(?![[:punct:]]|[[:alnum:]])"                 
# [5] "^(\\bimagination and bruises from playing)(?![[:punct:]]|[[:alnum:]])"     
# [6] "^(\\bwith imagination and bruises from playing)(?![[:punct:]]|[[:alnum:]])"

x<-"in almost all of Adam Sandler's"
pattern <-prepare(x)
pattern
# [1] "^(\\bsandler\\'s)(?![[:punct:]]|[[:alnum:]])"                      
# [2] "^(\\badam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"                 
# [3] "^(\\bof adam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"              
# [4] "^(\\ball of adam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"          
# [5] "^(\\balmost all of adam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"   
# [6] "^(\\bin almost all of adam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"

x<-"in almost all of Adam Sandler's"
pattern <-prepare(x)
pattern
# [1] "^(\\bsandler\\'s)(?![[:punct:]]|[[:alnum:]])"                      
# [2] "^(\\badam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"                 
# [3] "^(\\bof adam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"              
# [4] "^(\\ball of adam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"          
# [5] "^(\\balmost all of adam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"   
# [6] "^(\\bin almost all of adam sandler\\'s)(?![[:punct:]]|[[:alnum:]])"

#########


# sanity check on datfiles (full and sampled)
filename<-paste(datasetsdir,"FullData.RData",sep="/");load(file=filename)

fetchf<-grep(pattern[6],EN_data,value=FALSE,perl=FALSE,useBytes=TRUE);fx<-EN_data[fetchf]
fetchf<-grep(pattern[5],EN_data,value=FALSE,perl=FALSE,useBytes=TRUE);fx<-EN_data[fetchf]
fetchf<-grep(pattern[4],EN_data,value=FALSE,perl=FALSE,useBytes=TRUE);fx<-EN_data[fetchf]
fetchf<-grep(pattern[3],EN_data,value=FALSE,perl=FALSE,useBytes=TRUE);fx<-EN_data[fetchf]
fetchf<-grep(pattern[2],EN_data,value=FALSE,perl=FALSE,useBytes=TRUE);fx<-EN_data[fetchf]
fetchf<-grep(pattern[1],EN_data,value=FALSE,perl=FALSE,useBytes=TRUE);fx<-EN_data[fetchf]

filename<-paste(datasetsdir,"SampleData.RData",sep="/");load(file=filename)
fetchs<-grep(pattern[6],EN_datas,value=FALSE,perl=TRUE,useBytes=TRUE);fx<-EN_datas[fetchs]
fetchs<-grep(pattern[5],EN_datas,value=FALSE,perl=TRUE,useBytes=TRUE);fx<-EN_datas[fetchs]
fetchs<-grep(pattern[4],EN_datas,value=FALSE,perl=TRUE,useBytes=TRUE);fx<-EN_datas[fetchs]
fetchs<-grep(pattern[3],EN_datas,value=FALSE,perl=TRUE,useBytes=TRUE);fx<-EN_datas[fetchs]
fetchs<-grep(pattern[2],EN_datas,value=FALSE,perl=TRUE,useBytes=TRUE);fx<-EN_datas[fetchs]
fetchs<-grep(pattern[1],EN_datas,value=FALSE,perl=TRUE,useBytes=TRUE);fx<-EN_datas[fetchs]
