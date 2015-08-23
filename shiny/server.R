# server.R
library(shiny)
shinyServer ( function(input, output,session) { 
  output$oid1 <- renderText(paste0(tags$strong("You typed: "), tags$i(input$id1)))
  # the output from predictors
  observe ( {
       
       validate(
             need((grep("((\\b|[[:punct:]]|\\') +)$",input$id1,ignore.case=TRUE,perl=TRUE,value=FALSE,useBytes=TRUE) > 0),
                  message = "waiting for space or punctuation...",output$oid1 )
       )
       mstr <- input$id1
       if(exists("Gram")){
        SGT <- SGT_Skipgram_Predictor(input$id1)
        GOOGLE <- Google_Suggest_Predictor(input$id1)
        output$my.prediction <- renderUI(selectInput("nwordInput", "",choices=SGT, selectize=FALSE, multiple=TRUE, size=5))
        output$google.prediction <- renderUI(selectInput("gwordInput", "",choices=GOOGLE, selectize=FALSE, multiple=TRUE, size=5))
        output$my.newwords <- renderUI(selectInput("mwordInput","",choices=newword, selectize=FALSE, multiple=TRUE, size=5))
        # re benchmark...
        output$newPlot <- renderPlot({autoplot(microbenchmark(GOOGLE,SGT),times=1)}) 
        output$count <- renderUI({(paste("New Words Learned: ",length(newword)))})
        output$gramsize <- renderUI({(paste("Dictionary:",length(Gram[[1]]),"Words"))}) 
        output$session <- renderUI({(paste("Suggestion Number:",mysession))})
    }
  } )  
  
  observeEvent(input$nwordInput,{
    updateTextInput(session,inputId = "id1", value=paste0(input$id1,input$nwordInput," ") )
  })
  observeEvent(input$gwordInput,{
    updateTextInput(session,inputId = "id1", value=paste0(input$id1,input$gwordInput," ") )
  })
  observeEvent(input$mwordInput,{
    updateTextInput(session,inputId = "id1", value=paste0(input$id1,input$mwordInput," ") )
  })
} ) 
