# global.R
library(magrittr)
library(stringi)
library(ggplot2)
library(microbenchmark)
library(compiler)
# Gram limit is Quadrigram in this app
newword <-NULL
mysession <- 0L                     # session counter

if (!exists("Gram")) {Gram <-readRDS("GramX.RDS") }
SGT_Skipgram <- function(x) {
  #' initialize
  #' @param x the input string
  #' @return a pattern list of regex encoded tokens in x 
  #' @examples 
  #' initialize(mystring)   initializes, tokenizes and Regex prepares a patern list for grep 
  #' retain x, from a copy in y build a Regex pattern list 
mysession <<- mysession +1L
  stopifnot(is.character(x))
  if (nchar(x) == 0) {return(x)} # empty request!
  #' To avoid :: calls 
  stri_split_boundaries <- stringi::stri_split_boundaries 
  stri_join <- stringi::stri_join  
  stri_flatten <- stringi::stri_flatten
  options <- stringi::stri_opts_brkiter( type="word", skip_word_none = TRUE, skip_word_number = FALSE )
  # clean input
  x %<>%  gsub("[^A-Za-z\' ]","",.) %>%  tolower %>%
    stri_split_boundaries(opts_brkiter=options) %>% unlist
  w <- x  # keep a backup
  x %<>% tail(3) # only use last 3 words!
  x %>% length -> len
  # Split into word tokens
  x <- sapply (1L:len, function(i) {stri_flatten(tail(x,i),collapse=" ")})
  x %<>%  gsub("\\'","\\\\\\'",.,ignore.case=TRUE,perl=TRUE,useBytes=TRUE) 
  p <- x <- sapply (1:(len+1L), function(i) {
    ifelse (i > 1L, paste0("^(\\b",x[[i-1L]],")(?=[[:space:]])(.*)"),paste0("^(\\b",x[[1]],")(?![[:punct:]]|[[:alnum:]])") ) } ) 
  #' natural language processor to suggest continuation word
  #' The regex encoded list in y is processed to return a list of suggestion strings x in decreasing scores
  #' saved pattern list in p
  # use default 'Stupid Backoff' factor
  s <- 0.4
  n <- length(p)
  z <- list()
  z <- lapply (1:n, function(i) { .subset2(Gram,i)[grep(p[[i]],names(Gram[[i]]),ignore.case=TRUE,perl=TRUE,value=FALSE,useBytes=TRUE)] } ) 
  z <- Filter(length,z) 
  ## Predict Algorithm starts here
  n <- length(z)  # update length
  # check to see if we learned a new word today! if so update the grams for this session only if new words learned is TRUE
  if (n == 0) { 
  newword <<-c(tail(w,1),newword)   # always put last on top!
          Gram[[1]]<<-c(Gram[[1]],setNames(1L,tail(w,1)))
          if (length(w) > 1) {Gram[[2]]<<-c(Gram[[2]],setNames(1L,paste(w[(length(w)-1L)],w[(length(w))],sep=" ")))}
          if (length(w) > 2) {Gram[[3]]<<-c(Gram[[3]],setNames(1L,paste(w[(length(w)-2L)],w[(length(w)-1L)],w[(length(w))],sep=" ")))}
          if (length(w) > 3) {Gram[[4]]<<-c(Gram[[4]],setNames(1L,paste(w[(length(w)-3L)],w[(length(w)-2L)],w[(length(w)-1L)],w[(length(w))],sep=" ")))}
  }
  if (n > 1) {
          z <- sapply(1:n, function(i) { if (i > 1) { s^(n-i)*z[[i]]/sum(z[[i]]) } else { s^(n-i)*z[[i]]/sum(Gram[[i]])} } ) 
          z <- sapply(1:n, function(i) { head(z[[i]],5) })         # only work on top of the lists
          z <- sapply(1:n, function(i) { setNames (z[[i]], sapply(1:length(z[[i]]), function(j) { 
                  sub(p[[i]],"\\2",names(z[[i]][j]),ignore.case=TRUE,perl=TRUE,useBytes=TRUE) } ) ) } )
          z <- sapply(1:n, function(i) { z[[i]] <- setNames (z[[i]], sapply(1:length(z[[i]]), function(j) {
                  gsub("(^ )","",names(z[[i]][j]),ignore.case=TRUE,perl=TRUE,useBytes=TRUE) } ) ) } ) }  
  z[[1]] <- head(Gram[[1]],5)/sum(Gram[[1]]) 
  z %>% unlist %>% sort %>% rev %>% names %>% unique %>% head(5) %>% as.vector -> x
  # just in case we could not match anything, we always backoff with the 5 most common unigrams!
  x
}            
SGT_Skipgram_Predictor <- cmpfun (SGT_Skipgram, options = list(optimize=3))

Google_Suggest <- function(x) {
  #' retrieve suggestion from suggestqueries.google.com
  #' @param x the input string (must be char and of length >0)
  #' @return a list of suggestion strings in x in decreasing score rank 
  #' @examples 
  #' google_predictor('you hadn't time to take a') retruns list x: [1] "class"   "bath"    "break"   "nap"     "picture"
  ## Prepare x and return pattern list in x
  stopifnot (is.character(x)) 
  if (nchar(x) == 0) {return(x)} # empty request!
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
  x
} 
Google_Suggest_Predictor<-cmpfun (Google_Suggest,options = list(optimize=3)) 
