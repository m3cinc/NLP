library(rbenchmark)
library(microbenchmark)
library(proftools)
library(compiler)


skip_ngram2C <- function(x, n = 1L, skip = 0L, cflag = 0L) { 
        #' skip_ngram2 
        #' @param x character 
        #' @param n integer 
        #' @param skip integer
        #' @param cflag integer , 1 if compressed string x is input
        #' @return n-gram in x 
        #' @examples 
        #' skip_ngram(as.character(citation())) produces a unigram and compresses it
        #' skip_ngram(string,.,3) produces a trigram of the string
        #' skip_ngram(string,2,3) produces a skip3-bigram 
        library(magrittr)
        library(stringi)
        stopifnot(is.character(x), is.numeric(n), is.finite(n), n > 0, is.numeric(skip), is.finite(skip), skip >= 0 ,cflag >=0)
        if (cflag > 0) {x %>% memDecompress("g",asChar=TRUE) %>% strsplit(.,"\n") %>% unlist -> x}      # decompress
        if (n == 1 && skip > 0) { skip<-0 } # ignore the skip
        # To avoid :: calls 
        stri_split_boundaries <- stringi::stri_split_boundaries 
        stri_join <- stringi::stri_join  
        options <- stringi::stri_opts_brkiter( type="word", skip_word_none = TRUE, skip_word_number = FALSE )
        tokens <- unlist(stri_split_boundaries(x, opts_brkiter=options)) # Split into word tokens
        len <- length(tokens)
        # if we didn't detect any words or number of tokens is less than n return empty vector 
        if (all(is.na(tokens)) || len < n + skip) { character(0)} else { m <- n-1 ;
        if ( skip == 0) { sapply(1:(len-m), function(i) stri_join(tokens[i:(i+m)], collapse = " "))} else {
                          sapply(1:(len-m-skip), function(i) stri_join(c(tokens[i],tokens[(i+skip+1):(i+skip+m)]), collapse = " ")) } }
} 

skip_ngram2C_Level3<-cmpfun(skip_ngram2C,options = list(optimize=3))

#
# A number of test cases follow to demonstrate

compressor <- function(x) {
        stopifnot(is.character(x))
        x %>% memCompress("g") %>% memDecompress("g",asChar=TRUE) %>% strsplit(.,"\n") %>%
                extract2(1) %>% identical(x) -> wflag
        if (wflag==TRUE) memCompress(x,"g") else x
}


test <- c("This scope of this report is to perform exploratory data analysis on the Coursera Swiftkey data ",
 "set. The Swiftkey data set contains English, German, Russian and Finnish language data sets for",
 " exploration. This milestone report focuses on exploring English language datasets within the Swiftkey ",
 "data set. This data set contains en_US.news.txt, en_US.blogs.txt, and en_US.twitter.txt files.")
#
Rprof("Rprof-mem.out", memory.profiling=TRUE)
x<-skip_ngram2C(test,n=4,skip=0)
Rprof(NULL)
summaryRprof("Rprof-mem.out", memory="both") 
p <- readProfileData(filename="Rprof-mem.out")
plotProfileCallGraph(p, style=google.style, score="total")
printProfileCallGraph(p, file = stdout(), percent = TRUE)
gcinfo(TRUE)
x<-skip_ngram2C(test,n=2,skip=0)
gcinfo(FALSE)
library(microbenchmark)
microbenchmark(x<-skip_ngram2(test,n=4,skip=0))
library(compiler)
skip_ngram2_Level0<-cmpfun(skip_ngram2,options = list(optimize=0))
skip_ngram2_Level1<-cmpfun(skip_ngram2,options = list(optimize=1))
skip_ngram2_Level2<-cmpfun(skip_ngram2,options = list(optimize=2))
skip_ngram2C_Level3<-cmpfun(skip_ngram2C,options = list(optimize=3))

library(microbenchmark)
x <- runif(100)
bench <- microbenchmark (skip_ngram2(test,n=4,skip=0),
                         skip_ngram2_Level0(test,n=4,skip=0),
                         skip_ngram2_Level1(test,n=4,skip=0),
                         skip_ngram2_Level2(test,n=4,skip=0),
                         skip_ngram2_Level3(test,n=4,skip=0))
bench
Rprof("Rprof-mem.out", memory.profiling=TRUE)
x<-skip_ngram2C_Level3(test,n=2,skip=0)
Rprof(NULL)
summaryRprof("Rprof-mem.out", memory="both") 
p <- readProfileData(filename="Rprof-mem.out")
plotProfileCallGraph(p, style=google.style, score="total")
printProfileCallGraph(p, file = stdout(), percent = TRUE)
gcinfo(TRUE)
x<-skip_ngram2(test,n=4,skip=0)
gcinfo(FALSE)


# [1] "This scope of this"                                   "scope of this report"                                
# [3] "of this report is"                                    "this report is to"                                   
# [5] "report is to perform"                                 "is to perform explatory"                             
# ...
# [47] "This data set contains"                               "data set contains en_US.news.txt"                    
# [49] "set contains en_US.news.txt en_US.blogs.txt"          "contains en_US.news.txt en_US.blogs.txt and"         
# [51] "en_US.news.txt en_US.blogs.txt and en_US.twitter.txt" "en_US.blogs.txt and en_US.twitter.txt files"
#
# test1 <- skip_ngram(test,2,4);test1
# [1] "This of"                           "scope this"                        "of report"                        
# [4] "this is"                           "report to"                         "is perform"                       
# [7] "to explatory"                      "perform data"                      "explatory analysis"               
# ...
# [46] "set data"                          "This set"                          "data contains"                    
# [49] "set en_US.news.txt"                "contains en_US.blogs.txt"          "en_US.news.txt and"               
# [52] "en_US.blogs.txt en_US.twitter.txt" "and files"                        
#
# test2 <- skip_ngram(test,4,2);test2
# [1] "This this report is"                          "scope report is to"                          
# [3] "of is to perform"                             "this to perform explatory"                   
# [5] "report perform explatory data"                "is explatory data analysis"                  
# [7] "to data analysis on"                          "perform analysis on the"                     
# ...
# [47] "This contains en_US.news.txt en_US.blogs.txt" "data en_US.news.txt en_US.blogs.txt and"     
# [49] "set en_US.blogs.txt and en_US.twitter.txt"    "contains and en_US.twitter.txt files"        
