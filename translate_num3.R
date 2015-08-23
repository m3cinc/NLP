#' translate3_num
#' @param x character
#' @param spell_number a named int list of spelled numbers [0-1000] 
#' @examples 
#' translate3_num( as.character( citation()), spell_number) substitutes number with spelled value
#' nskipgram(string) 

translate_numC <- cmpfun( translate_num <- function ( x, spell_number ) { 
        stopifnot( is.character( x ))
        #' To avoid :: calls 
        stri_split_boundaries <- stringi::stri_split_boundaries 
        stri_join <- stringi::stri_join  
        stri_flatten <- stringi::stri_flatten
        options <- stringi::stri_opts_brkiter( type = "word", skip_word_none = TRUE, skip_word_number = FALSE )
        if ( length( grep( "([#0-9])", x, perl = TRUE, fixed = FALSE, value = FALSE ) ) == 0 ) { return(x) } 
        # so we have possibly # or digits
        x %<>% gsub("(#)(\\d+)","\\1 \\2",.,perl = TRUE)  # a combination like "#12" translates into "# 12"
        x %<>% gsub("# "," number ",.)                    # a combination like "# " translates into "number "
        x %<>% gsub("(\\d+)","\\1 ",.,perl = TRUE)        # a combination like "2go" translated to " 2 go "
        x %<>% gsub("(#\\S+)","#hashtag ",.,perl = TRUE)  # a hashtag converted to "#hashtag "
        if ( length( grep( "([0-9])",x , perl = TRUE, fixed = FALSE, value = FALSE ) ) == 0 ) { 
                # just flatten before return
                x <- stri_flatten(x, collapse = " ") 
                return(x) }
        tokens <- unlist( stri_split_boundaries(x, opts_brkiter=options)) # Split into word tokens
        # If we didn't detect any words or number of tokens is less than n return empty vector 
        if (length(tokens) == 1L) { x <- tokens[length(tokens)] ; return(x) }
        change_tokens <- grep( "([0-9]){1,3}", tokens, perl = TRUE, value = FALSE )
        if (length(change_tokens) == 0L) {return(x)}
        tokens[change_tokens] <- sapply(change_tokens, function(i) { 
                tokens[i] <- ifelse (as.integer(grep( "([0-9]){1,3}", tokens[i], perl = TRUE, value = TRUE ) ) %in% 0L:999L,
                                      gsub( "([0-9]){1,3}", as.character( names( spell_number[1L+as.integer(tokens[i])])),tokens[i],perl=TRUE),tokens[i])
        } )  
        x <- stri_flatten( tokens, collapse = " ") 
} , options = list( optimize = 3) )
  
#
# A number of test cases follow to demonstrate
#
#test <- c("This scope of this report is to perform explatory data analysis on #101 pages of the Coursera Swiftkey data ",
#          "set. The Swiftkey'data set contains English, German, Russian and Finnish language data sets for",
#          " exploration. This milestone report # 213 focuses on exploring English language datasets within the Swiftkey ",
#          "data set. This data set contains en_US.news.txt, en_US.blogs.txt, and en_US.twitter.txt files.")
