# ui.R
library(shiny)
shinyUI (fluidPage (
        titlePanel("SGT Skip-Gram Turing NLP Predictor"),
        sidebarLayout(
          sidebarPanel(
                h5('About Me ...'),
                h5('My task is to predict 5 next words to continue a sentence You type, ranked-ordered from most- to least-likely.
                    I will also show you what trained cousin GoogleSuggest(TM) predicts!'),
                h5('As a Shiny App, I leverage:'),
                h5('1.  ',tags$a(href="https://class.coursera.org/nlp/lecture/128", target="_blank", "Quadrigram Vocabulary"),' to maximize coverage,'),
                h5('2.  ',tags$a(href="https://faculty.cs.byu.edu/~ringger/CS479/papers/Gale-SimpleGoodTuring.pdf",target="_blank","Simple Good Turing"),' smoothing and pruning,'),
                h5('3.  ',tags$a(href="http://homepages.inf.ed.ac.uk/ballison/pdf/lrec_skipgrams.pdf",target="_blank","Skip-Gram"),' prediction range extender,'),
                h5('4.  ',span ("Learning Mode", style = "color:blue"),' to demonstrate in-session capabilities,'),
                h5('5.  ',tags$a(href="https://regex101.com/",target="_blank","Regex-based"),' filtering, RRRepetition and noise suppressor,'),
                h5('6.  ',tags$a(href="https://en.wikipedia.org/wiki/Katz%27s_back-off_model", target="_blank","Proven Katz' Backoff"),' NLP model algorithm,'),
                h5('7.  Trained cousin ',tags$a(href="https://support.google.com/websearch/answer/106230",target="_blank","GoogleSuggest(TM)"),' suggestion,'),
                h5('8.  ',tags$a(href="http://adv-r.had.co.nz/Profiling.html", target="_blank","Compiler and Profiler"),' performance enhancers,'),
                h5('9.  ',span ("Dictionary", style = "color:blue"),' to fit mobile and fixed application useage,'),
                h5('... And much more, as you are bound to discover...'),
                h5('As a Shiny App, I display ', span (" Benchmarked Performance Statistics", style = "color:blue"),' on the fly'),
                h5("I was trained on a randomized 80% subset of the 550 MBytes Swift data, representing all the filtered 
                    twitter, blogs and news my builder's hardware could sustain, but will learn more vocabulary as I go..."),
                h5('I was developed for the Summer 2015 Edition of the',tags$a(href="https://www.coursera.org/specialization/jhudatascience/1",target="_blank",
                        "Data Science Coursera track"),'Capstone project, by a Data Science novice who had a lot of fun learning and practicing on a production laptop!'),
                h5(''),
                h5('I am packed with Innovation, Technology, Benchmark and Performance Features to empower You!'),
                h5('I am committed to, and will continue to improve! Have Fun and Enjoy!'),
                fluidRow(
                        column(3, uiOutput("gramsize")),
                        column(3, uiOutput("count")),
                        column(3, uiOutput("session"))
                )        
                
          ), 
          mainPanel(
                conditionalPanel("mysession.length > 0"),
                
                fluidRow(
                    column(12,
                        tags$h5("You type your starting sentence. I will predict next words and display in ordered lists below as soon as I detect a space in your input ..."),
                        tags$h5("You may then proceed by typing or picking from any of the predicted lists... Let's have Fun!"),
 #                       div(style="display:inline-block", submitButton("Next Word!"),width=4),
                        div(style="display:inline-block",textInput("id1", "",value="I will do this once, so you see what I ")), 
                        tags$head(tags$style(type="text/css", "#id1 {width: 1020px}"))
                  )
                ),
                fluidRow(
                    column(12,    
                    uiOutput("oid1")
                    )
                ),
                fluidRow(
                  column(4,
                    h5('SGT predicts the next word will be ...'),
                    uiOutput("my.prediction")
                  ),
                  column(4,
                    h5('Cousin GoogleSuggest(TM) predicts ...'),
                    uiOutput("google.prediction")
                  ),
                  column(4,
                    h5('SGT New Learned Words list is ...'),
                    uiOutput("my.newwords")
                  )
                ),
 #               tags$hr(),
 #               div(style="display:inline-block", submitButton("Pick!"),width=4),
 #               div(style="display:inline-block",tags$h4("a continuing word in either lists above or keep typing...")),
 #               tags$hr(),
 #               h5('While you type in another sentence, I will tally the benchmark to keep track of my performance!'),
                plotOutput("newPlot")
         )
        )
)) 
