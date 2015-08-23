#library(rbenchmark)
#library(microbenchmark)
#library(proftools)
library(compiler)
#library(ggplot2)
skip_ngram3C <- cmpfun( skip_ngram3 <- function( x, badwords, n = 1L, skip = 0L, cflag = 0L) { 
        #' skip_ngram2 
        #' @param x character 
        #' @param n integer 
        #' @param skip integer
        #' @param cflag integer, 1L if compressed string x is input
        #' @return n-gram in x 
        #' @examples 
        #' skip_ngram3C(as.character(citation())) produces a unigram of citation in list x
        #' skip_ngram3C(string,.,3) produces trigram of string in list x
        #' skip_ngram3C(string,2,3) produces skip3-bigram in list x corresponding to [word XXX word]
        stopifnot( is.integer( n ), is.finite( n ), n > 0, is.integer( skip ), is.finite( skip ), skip >-1L, is.integer( cflag ), is.finite( cflag ), cflag >-1L)
#        require(magrittr,stringi)
        if ( cflag > 0 ) { x %>% memDecompress( "g", asChar = TRUE ) %>% strsplit( "\n" ) %>% unlist -> x }      # decompress
        if ( n == 1 && skip > 0 ) { skip <-0 } # ignore the skip
        # To avoid :: calls 
        stri_split_boundaries <- stringi::stri_split_boundaries 
        stri_join <- stringi::stri_join  
        options <- stringi::stri_opts_brkiter( type = "word", skip_word_none = TRUE, skip_word_number = FALSE )
        tokens <- unlist( stri_split_boundaries( x, opts_brkiter = options ) ) # Split into word tokens
        tokens <- tokens[!tokens %in% badwords]   # eliminate undesirable vocabulary
        # if we didn't detect any words or number of tokens is less than n return empty vector 
        if ( all( is.na( tokens )) || length( tokens ) < n + skip ) { character(0) } else { 
        if ( skip == 0) { sapply(1:( length( tokens ) - n - 1L),       function(i) stri_join( tokens[i:(i + n - 1L)], collapse = " "))} else {
                          sapply(1:( length( tokens ) - n - skip -1L ), function(i) stri_join( c( tokens[i], tokens[(i + skip + 1L):(i + skip + n - 1L)]), collapse = " ")) } }
} , options = list( optimize = 3) )

#
# A number of test cases follow to demonstrate

#compressor_check <- function( x ) {
#        stopifnot( is.character( x ))
#        require( magrittr )
#        x %>% memCompress( "g" ) %>% memDecompress( "g", asChar = TRUE ) %>% strsplit( "\n" ) %>% magrittr::extract2(1) %>% identical(x) -> wflag
#        if ( wflag == TRUE ) { memCompress( x, "g" ) } else x
#}
#
#
# test <- c( "This scope of this report is to perform exploratory data analysis on the Coursera Swiftkey data ",
# "set. The Swiftkey data set contains English, German, Russian and Finnish language data sets for",
# " exploration. This milestone report focuses on exploring English language datasets within the Swiftkey ",
# "data set. This data set contains en_US.news.txt, en_US.blogs.txt, and en_US.twitter.txt files." )
#
#Rprof( "Rprof-mem.out", memory.profiling = TRUE)
#x <- skip_ngram3( test, n = 4L , skip = 0L )
#Rprof( NULL)
#summaryRprof( "Rprof-mem.out", memory = "both" ) 
#p <- readProfileData( filename = "Rprof-mem.out" )
#plotProfileCallGraph( p, style = google.style, score = "total")
#printProfileCallGraph( p, file = stdout(), percent = TRUE)
#gcinfo( TRUE)
#x <- skip_ngram3(test, n = 2L, skip = 0L)
#gcinfo( FALSE)
#microbenchmark( x <- skip_ngram3( test, n = 4L, skip = 0L) )
#skip_ngram3_Level0 <- cmpfun( skip_ngram3, options = list( optimize = 0))
#skip_ngram3_Level1 <- cmpfun( skip_ngram3, options = list( optimize = 1))
#skip_ngram3_Level2 <- cmpfun( skip_ngram3, options = list( optimize = 2))
#skip_ngram3_Level3 <- cmpfun( skip_ngram3, options = list( optimize = 3))
#
#x <- runif(100)
#Rprof( "Rprof-mem.out", memory.profiling = TRUE)
#bench <- microbenchmark ( skip_ngram3( test, n = 4L, skip = 0L),
#                          skip_ngram3_Level0( test, n = 4L, skip = 0L),
#                          skip_ngram3_Level1( test, n = 4L, skip = 0L),
#                          skip_ngram3_Level2( test, n = 4L, skip = 0L),
#                          skip_ngram3_Level3( test, n = 4L, skip = 0L))
#autoplot(bench)
#Rprof( "Rprof-mem.out", memory.profiling = TRUE)
#x <- skip_ngram3C( test, n = 2L, skip = 0L)
#Rprof( NULL)
#summaryRprof( "Rprof-mem.out", memory = "both") 
#p <- readProfileData( filename = "Rprof-mem.out")
#plotProfileCallGraph( p, style = google.style, score = "total")
#printProfileCallGraph( p, file = stdout(), percent = TRUE)
#gcinfo( TRUE)
#x <- skip_ngram3C( test, n = 4L, skip = 0L)
#gcinfo( FALSE)
#autoplot(microbenchmark(skip_ngram3(test, n = 2L, skip = 0L),skip_ngram3C(test, n = 2L, skip = 0L),times=50))
#
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
