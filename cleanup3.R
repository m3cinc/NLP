#' cleanup3.R
#' @param x character
#' @param badwords a list of bad words characters
#' @param spell_number a maned int to spell numbers 0-1000
#' @examples 
#' cleanup(as.character(citation()),badwords,spell_number) cleans and substitute citation()
#' cleanup(string,badwords,spell_number) 

cleanup3C <- cmpfun ( cleanup3 <- function (x, spell_number) { 
        stopifnot(is.character(x))
        require(magrittr,stingi)
        x %>% 
                gsub("[^A-Za-z\\' 0-9$+=@#%^&_<>*.-]","",.) %>%
                gsub("[\\w\\-][\\w\\-\\.]+@[\\w\\-][\\w\\-\\.]+[a-zA-Z]{1,4}","email-address",.,perl=TRUE) %>%
                gsub("(?:(?:https?|ftp|file):\\/\\/|www\\.|ftp\\.)(?:\\([-A-Z0-9+&@#\\\\/%=~_|$?!:;,.]*\\)|[-A-Z0-9+&@#\\\\/%=~_|$?!:;,.])*(?:\\([-A-Z0-9+&@#\\\\/%=~_|$?!:;,.]*\\)|[A-Z0-9+&@#\\\\/%=~_|$](?x))","webaddress",.,perl=TRUE,ignore.case=TRUE) %>%
                gsub("( +([Uu] +[Rr]( *|\\.))+)"," you are ",.,perl=TRUE) %>%
                gsub("( +([Uu][Rr]( |\\.))+)"," your\\3",.,perl=TRUE) %>%
                gsub("( +([Rr][Uu]( |\\.))+)"," are you\\3",.,perl=TRUE) %>%
                gsub("( +([Ll][Uu]( |\\.))+)"," love you\\3",.,perl=TRUE) %>%
                gsub("( +([Rr][Tt]( |\\.))+)","  return\\3",.,perl=TRUE) %>%
                gsub("( +([Rr]( |\\.))+)"," are\\3",.,perl=TRUE) %>%
                gsub("( +([Uu]( |\\.))+)"," you\\3",.,perl=TRUE) %>%
                gsub("( +([bb]( |\\.))+)"," be\\3",.,perl=TRUE) %>%
                gsub("( +([Yy]( |\\.))+)"," why\\3",.,perl=TRUE) %>%
                gsub("((\\$)(\\d+((\\.|\\,)\\d+))+)","\\3 dollars",.,perl=TRUE) %>% # $ in dollars
                gsub("(((\\$\\.)(\\d+))+)","\\4 cents",.,perl=TRUE) %>%             # fractional $ in cents
                gsub("((\\d+)(\\% )+)","\\2 percent",.,perl=TRUE) %>%               # percent
                gsub("(\\&|\\+)( *[A-Za-z]+)","and \\2",.,perl=TRUE) %>%            # (&|+)[[:alpha:]] in and [[:alpha:]]
                gsub("(\\&|\\+)( *[A-Za-z]+)","plus \\2",.,perl=TRUE) %>%           # (&|+)[[:digits:]] in plus [[:digits:]]
                gsub("(\\-)( *[0-9]+)","minus \\2 ",.,perl=TRUE) %>%                # -[[:digits]] in minus [[:digits:]]
                gsub("(\\-)( *[A-Za-z]+)","less \\2 ",.,perl=TRUE) %>%              # -[[:alpha]] in less [[:alpha:]]
                gsub("(\\=)( *[A-Za-z]+)","equal ",.,perl=TRUE) %>%                 # = in equal
                gsub("(\\<)( *[A-Za-z]+)","less than ",.,perl=TRUE) %>%             # < in less than
                gsub("(\\>)( *[A-Za-z]+)","greater than ",.,perl=TRUE) %>%          # > in greater than
                stri_split_boundaries(.,n=-1L,type="sentence" ) %>%  unlist %>% stri_split_lines(omit_empty=TRUE) %>% stri_trim_both(pattern="\\P{Wspace}") %>% gsub("\\s+"," ",.,perl=TRUE) %>% unlist -> x
        x[x != "" & !is.na(x)] ->x
        x %>% grep("( +[#0-9]){1,3}",.,perl=TRUE,fixed=FALSE,value=FALSE) -> subline
        if (length(subline) > 0) { x[subline]<-pbsapply(subline, function(i) {translate_numC(x[(i)], spell_number)}) }
        x %>% gsub("[^A-Za-z \\'\\-]|[0-9]","",.,perl=TRUE) %>%  stri_trim_both(pattern="\\P{Wspace}") %>% gsub("\\s+"," ",.,perl=TRUE) %>% unlist %>% tolower %>% gsub("( +([^ai])( +|\\.))|(^[^ai])( +|\\.)|( +[^ai]$)( *|\\.)"," ",.,perl=TRUE) -> x
        x[x != "" & !is.na(x)]-> x
}, options = list( optimize = 3) )
