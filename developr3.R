# developr3.R
options(mc.cores=4)
memory.limit(size=28000)
setwd("V:/") 
Sys.info()[1:5]                         # system info but exclude login and user info
userdir<-getwd()                        # user-defined startup Capstone directory
library(magrittr)
library(pbapply)
library(stringi)
library(compiler)
local <- TRUE                           # to avoid data redownload
datadir <- "./data"
datasetsdir <- paste0(datadir,"/datasets")
userdatadir <- paste0(datadir,"/Final/en_US")
gramdir <- paste0(datadir,"/gram")
elimdir <- paste0(datadir,"/elim")

if (!file.exists(datadir)) { dir.create(datadir) }  #  data will reside in datadir
if (!file.exists(datasetsdir)) { dir.create(datasetsdir) }  #  datasets will reside in datasetsdir
if (!file.exists(gramdir)) { dir.create(gramdir) }  #  gram will reside in datadir
if (!file.exists(elimdir)) { dir.create(elimdir) }  #  gram will reside in datadir
###############################################
#
# read raw data
#
###############################################
US.filenames <- list.files(path=userdatadir,full.names=TRUE)
US.blogs <- US.filenames[[1]]
US.news <- US.filenames[[2]]
US.twitter <- US.filenames[[3]]
con <- file( US.blogs, open = "rb", encoding = "UTF-8" ) ; blogs_raw <- readLines( con, encoding = "UTF-8" ) ; close( con )
con <- file( US.news, open = "rb", encoding = "UTF-8" ) ; news_raw <- readLines( con, encoding = "UTF-8" ) ; close( con)
con <- file( US.twitter, open = "rb", encoding = "UTF-8" ) ; twitter_raw <- readLines( con, encoding = "UTF-8" ) ; close( con )
rm(con, US.blogs, US.filenames, US.news, US.twitter)
###############################################
#
# step 1. process into ASCII
#
###############################################
blogs_raw   %<>% iconv( "latin1", "ASCII", sub = "")
news_raw    %<>% iconv( "latin1", "ASCII", sub = "")
twitter_raw %<>% iconv( "latin1", "ASCII", sub = "")
###############################################
#
# step 2. substitute numbers, and other symbols (+-$%#@&...)
#
###############################################
filename <- paste( datadir, "en_US.badwords.txt", sep = "/" ); con <- file( filename, open = "rb" ) ; badwords <- readLines(con, encoding = "UTF-8") ; close(con)
EN_data <- c( blogs_raw , news_raw , twitter_raw )
saveRDS( EN_data, file = "EN_data_US.RDS" )
rm ( blogs_raw , news_raw , twitter_raw ,con)
filename <- paste0(userdir,"spell3.R") ; source(filename)   # builds the spell_number named integer array [0L:1000L]
filename <- paste0( userdir, "translate_num3.R" ) ; source ( filename )
filename <- paste0( userdir, "cleanup3.R" ) ; source( filename )
EN_data %<>% cleanup3C( spell_number )     # clean and spell numbers [0L:1000L], and some other filters (%,$,#,@, etc...)
###############################################
#
# step 3. randomize sampling into train(50%), hold(25%) and test(25%) sets 
#
###############################################
rbsampler <- function(x,p) { return(x[as.logical( rbinom( length( x ), 1, p))]) }
set.seed(123123)                                        # set the seed to make your partition reproductible
train <- rbsampler(EN_data,0.5)                         # 50% of the sample size
hold <- rbsampler(EN_data[!EN_data %in% train] ,0.5)    # 50% of the remainder
test <- EN_data[!EN_data %in% train & !EN_data %in% hold]
saveRDS(train,"train.RDS") ; saveRDS(test,"test.RDS") ; saveRDS(EN_data,"EN_data.RDS")
EN_data <- train                                        # work on the train set
EN_data %<>% memCompress("g")                           # we need to save memory for making grams
rm( test, hold, train, con, spell_number, userdatadir )                          
###############################################
#
# step 4. prepare grams 
#
###############################################
filename<-paste0(userdir,"skip_ngram3.R") ; source(filename)
# make unigram & sanitize
EN_data %>% skip_ngram3C( badwords,1L,0L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram10 ; filename <- paste( gramdir, "S50Gram10.RData", sep = "/" ) ; save ( list = c( "S50Gram10" ), file = filename)
S50Gram10 %>% extract(.<2) %>% names -> Discard10 ; filename <- paste( elimdir, "Discard10.RDS", sep = "/" ) ; saveRDS( Discard10, file = filename) ; rm( Discard10 )
S50Gram10 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram10.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram10 )  # trimmed S50Gram10
# make bigram & sanitize
EN_data %>% skip_ngram3C( badwords,2L,0L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram20 ; filename <- paste( gramdir, "S50Gram20.RData", sep = "/" ) ; save ( list = c( "S50Gram20" ), file = filename)
S50Gram20 %>% extract(.<2) %>% names -> Discard20 ; filename <- paste( elimdir, "Discard20.RDS", sep = "/" ) ; saveRDS( Discard20, file = filename) ; rm( Discard20 )
S50Gram20 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram20.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram20 )  # trimmed S50Gram20
# make trigram & sanitize
EN_data %>% skip_ngram3C( badwords,3L,0L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram30 ; filename <- paste( gramdir, "S50Gram30.RData", sep = "/" ) ; save ( list = c( "S50Gram30" ), file = filename)
S50Gram30 %>% extract(.<2) %>% names -> Discard30 ; filename <- paste( elimdir, "Discard30.RDS", sep = "/" ) ; saveRDS( Discard30, file = filename) ; rm( Discard30 )
S50Gram30 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram30.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram30 )  # trimmed S50Gram30
# make quadrigram & sanitize
EN_data %>% skip_ngram3C( badwords,4L,0L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram40 ; filename <- paste( gramdir, "S50Gram40.RData", sep = "/" ) ; save ( list = c( "S50Gram40" ), file = filename)
S50Gram40 %>% extract(.<2) %>% names -> Discard40 ; filename <- paste( elimdir, "Discard40.RDS", sep = "/" ) ; saveRDS( Discard40, file = filename) ; rm( Discard40 )
S50Gram40 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram40.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram40 )  # trimmed S50Gram40
# make quintigram & sanitize
EN_data %>% skip_ngram3C( badwords,5L,0L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram50 ; filename <- paste( gramdir, "S50Gram50.RData", sep = "/" ) ; save ( list = c( "S50Gram50" ), file = filename)
S50Gram50 %>% extract(.<2) %>% names -> Discard50 ; filename <- paste( elimdir, "Discard50.RDS", sep = "/" ) ; saveRDS( Discard50, file = filename) ; rm( Discard50 )
S50Gram50 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram50.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram50 )  # trimmed S50Gram50
# make sextigram & sanitize
EN_data %>% skip_ngram3C( badwords,6L,0L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram60 ; filename <- paste( gramdir, "S50Gram60.RData", sep = "/" ) ; save ( list = c( "S50Gram60" ), file = filename)
S50Gram60 %>% extract(.<2) %>% names -> Discard60 ; filename <- paste( elimdir, "Discard60.RDS", sep = "/" ) ; saveRDS( Discard60, file = filename) ; rm( Discard60 )
S50Gram60 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram60.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram60 )  # trimmed S50Gram60
# make septigram & sanitize
EN_data %>% skip_ngram3C( badwords,7L,0L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram70 ; filename <- paste( gramdir, "S50Gram70.RData", sep = "/" ) ; save ( list = c( "S50Gram70" ), file = filename)
S50Gram70 %>% extract(.<2) %>% names -> Discard70 ; filename <- paste( elimdir, "Discard70.RDS", sep = "/" ) ; saveRDS( Discard70, file = filename) ; rm( Discard70 )
S50Gram70 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram70.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram70 )  # trimmed S50Gram70
###############################################
#
# step 5. prepare skip grams 
#
###############################################
# make bigrams & sanitize
EN_data %>% skip_ngram3C( badwords,2L,1L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram21 ; filename <- paste( gramdir, "S50Gram21.RData", sep = "/" ) ; save ( list = c( "S50Gram21" ), file = filename)
S50Gram21 %>% extract(.<2) %>% names -> Discard21  ; filename <- paste( elimdir, "Discard21.RDS", sep = "/" ) ; saveRDS( Discard21, file = filename) ; rm( Discard21 )
S50Gram21 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram21.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram21 )  # trimmed S50Gram21
EN_data %>% skip_ngram3C( badwords,2L,2L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram22 ; filename <- paste( gramdir, "S50Gram22.RData", sep = "/" ) ; save ( list = c( "S50Gram22" ), file = filename)
S50Gram22 %>% extract(.<2) %>% names -> Discard22  ; filename <- paste( elimdir, "Discard22.RDS", sep = "/" ) ; saveRDS( Discard22, file = filename) ; rm( Discard22 )
S50Gram22 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram22.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram22 )  # trimmed S50Gram22
EN_data %>% skip_ngram3C( badwords,2L,3L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram23 ; filename <- paste( gramdir, "S50Gram23.RData", sep = "/" ) ; save ( list = c( "S50Gram23" ), file = filename)
S50Gram23 %>% extract(.<2) %>% names -> Discard23  ; filename <- paste( elimdir, "Discard23.RDS", sep = "/" ) ; saveRDS( Discard23, file = filename) ; rm( Discard23 )
S50Gram23 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram23.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram23 )  # trimmed S50Gram23
EN_data %>% skip_ngram3C( badwords,2L,4L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram24 ; filename <- paste( gramdir, "S50Gram24.RData", sep = "/" ) ; save ( list = c( "S50Gram24" ), file = filename)
S50Gram24 %>% extract(.<2) %>% names -> Discard24  ; filename <- paste( elimdir, "Discard24.RDS", sep = "/" ) ; saveRDS( Discard24, file = filename) ; rm( Discard24 )
S50Gram24 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram24.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram24 )  # trimmed S50Gram24
# make trigrams & sanitize
EN_data %>% skip_ngram3C( badwords,3L,1L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram31 ; filename <- paste( gramdir, "S50Gram31.RData", sep = "/" ) ; save ( list = c( "S50Gram31" ), file = filename)
S50Gram31 %>% extract(.<2) %>% names -> Discard31 ; filename <- paste( elimdir, "Discard21.RDS", sep = "/" ) ; saveRDS( Discard31, file = filename) ; rm( Discard31 )
S50Gram31 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram31.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram31 )  # trimmed S50Gram31
EN_data %>% skip_ngram3C( badwords,3L,2L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram32 ; filename <- paste( gramdir, "S50Gram32.RData", sep = "/" ) ; save ( list = c( "S50Gram32" ), file = filename)
S50Gram32 %>% extract(.<2) %>% names -> Discard32 ; filename <- paste( elimdir, "Discard22.RDS", sep = "/" ) ; saveRDS( Discard32, file = filename) ; rm( Discard32 )
S50Gram32 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram32.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram32 )  # trimmed S50Gram32
EN_data %>% skip_ngram3C( badwords,3L,3L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram33 ; filename <- paste( gramdir, "S50Gram33.RData", sep = "/" ) ; save ( list = c( "S50Gram33" ), file = filename)
S50Gram33 %>% extract(.<2) %>% names -> Discard33 ; filename <- paste( elimdir, "Discard23.RDS", sep = "/" ) ; saveRDS( Discard33, file = filename) ; rm( Discard33 )
S50Gram33 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram33.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram33 )  # trimmed S50Gram33
# make quadrigrams & sanitize
EN_data %>% skip_ngram3C( badwords,4L,1L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram41 ; filename <- paste( gramdir, "S50Gram41.RData", sep = "/" ) ; save ( list = c( "S50Gram41" ), file = filename)
S50Gram41 %>% extract(.<2) %>% names -> Discard41 ; filename <- paste( elimdir, "Discard41.RDS", sep = "/" ) ; saveRDS( Discard41, file = filename) ; rm( Discard41 )
S50Gram41 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram41.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram41 )  # trimmed S50Gram41
EN_data %>% skip_ngram3C( badwords,4L,2L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram42 ; filename <- paste( gramdir, "S50Gram42.RData", sep = "/" ) ; save ( list = c( "S50Gram42" ), file = filename)
S50Gram42 %>% extract(.<2) %>% names -> Discard42 ; filename <- paste( elimdir, "Discard42.RDS", sep = "/" ) ; saveRDS( Discard42, file = filename) ; rm( Discard42 )
S50Gram42 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram42.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram42 )  # trimmed S50Gram42
# make quintigrams & sanitize
EN_data %>% skip_ngram3C( badwords,5L,1L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram51 ; filename <- paste( gramdir, "S50Gram51.RData", sep = "/" ) ; save ( list = c( "S50Gram51" ), file = filename)
S50Gram51 %>% extract(.<2) %>% names -> Discard51 ; filename <- paste( elimdir, "Discard51.RDS", sep = "/" ) ; saveRDS( Discard51, file = filename) ; rm( Discard51 )
S50Gram51 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram51.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram51 )  # trimmed S50Gram51

EN_data %>% skip_ngram3C( badwords,2L,5L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram25 ; filename <- paste( gramdir, "S50Gram25.RData", sep = "/" ) ; save ( list = c( "S50Gram25" ), file = filename)
S50Gram25 %>% extract(.<2) %>% names -> Discard25  ; filename <- paste( elimdir, "Discard25.RDS", sep = "/" ) ; saveRDS( Discard25, file = filename) ; rm( Discard25 )
S50Gram25 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram25.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram25 )  # trimmed S50Gram25
EN_data %>% skip_ngram3C( badwords,3L,4L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram34 ; filename <- paste( gramdir, "S50Gram34.RData", sep = "/" ) ; save ( list = c( "S50Gram34" ), file = filename)
S50Gram34 %>% extract(.<2) %>% names -> Discard34  ; filename <- paste( elimdir, "Discard34.RDS", sep = "/" ) ; saveRDS( Discard34, file = filename) ; rm( Discard34 )
S50Gram34 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram34.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram34 )  # trimmed S50Gram34
EN_data %>% skip_ngram3C( badwords,4L,3L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram43 ; filename <- paste( gramdir, "S50Gram43.RData", sep = "/" ) ; save ( list = c( "S50Gram43" ), file = filename)
S50Gram43 %>% extract(.<2) %>% names -> Discard43 ; filename <- paste( elimdir, "Discard43.RDS", sep = "/" ) ; saveRDS( Discard43, file = filename) ; rm( Discard43 )
S50Gram43 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram43.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram43 )  # trimmed S50Gram43
EN_data %>% skip_ngram3C( badwords,5L,2L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram52 ; filename <- paste( gramdir, "S50Gram52.RData", sep = "/" ) ; save ( list = c( "S50Gram52" ), file = filename)
S50Gram52 %>% extract(.<2) %>% names -> Discard52 ; filename <- paste( elimdir, "Discard52.RDS", sep = "/" ) ; saveRDS( Discard52, file = filename) ; rm( Discard52 )
S50Gram52 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram52.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram52 )  # trimmed S50Gram52
EN_data %>% skip_ngram3C( badwords,6L,1L,1L) %>% unlist %>% table %>% sort(decreasing=TRUE) -> S50Gram61 ; filename <- paste( gramdir, "S50Gram61.RData", sep = "/" ) ; save ( list = c( "S50Gram61" ), file = filename)
S50Gram61 %>% extract(.<2) %>% names -> Discard61 ; filename <- paste( elimdir, "Discard61.RDS", sep = "/" ) ; saveRDS( Discard61, file = filename) ; rm( Discard61 )
S50Gram61 %>% extract(.>1) -> Vocabulary ; filename <- paste( gramdir, "CS50Gram61.RData", sep = "/" ) ; save( list = c( "Vocabulary" ), file = filename ) ; rm( S50Gram61 )  # trimmed S50Gram61
