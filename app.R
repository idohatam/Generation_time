
library(shiny)
library(dplyr)
library(tidyr)
library(plotly)
library(shinycssloaders)

gnt3 <- function(g = sample(seq(20,1500, by=10),1), 
                 b = sample(seq(100,10000, by = 10),1)){
    #vector for generation time
    t <- (2:16)*g
    #vector for cell numlog phase
    cl <- rep(0, 15)
    
    cl_lag <- rep(0,100)
    
    #time lag
    
    t_lag <- seq(0,t[1], along.with = cl_lag)
    
    
    
    i=1
    
    for(i in seq_along(cl_lag)){
        if(i == 1){
            cl_lag[i] <- b
        }else{
            if(i > 1){
                cl_lag[i] <- cl_lag[i-1]+(cl_lag[i-1]*0.001)
            }
        }
    }
    
    i=1
    
    #loop to get cell numbers of log phase
    for (i in seq_along(cl)) {
        if(i==1){
            cl[i] <- cl_lag[length(cl_lag)]+(cl_lag[length(cl_lag)]*0.001)
        } else {
            
            if(i > 1){
                cl[i] <- cl[i-1]*2
            }
        }
        
    }
    
    #log transform clle numbers
    cl <- log10(cl)
    # linear regretion between time and cell number
    
    lnr <- lm(t~cl)
    
    #extrapulation
    cl2 <- seq(min(cl),max(cl), by = 0.0001)
    
    t2 <- ((cl2*lnr$coefficients[2])+lnr$coefficients[1])
    
    
    #get stat
    cl_stat<- rep(0,100)
    i=1
    for(i in seq_along(cl_stat)){
        if(i==1){
            cl_stat[i] <- cl2[length(cl2)]+(cl2[length(cl2)]*0.00005)
        } else{
            if(i > 1){
                cl_stat[i] <- cl_stat[i-1] + (cl_stat[i-1]*0.00005)
            }
        }
    }
    
    
    #get time stat
    t_stat <- seq(t2[length(t2)], t2[length(t2)]+(2*g), along.with = cl_stat)
    
    
    #get cell vector
    
    cl_tot <- c(log10(cl_lag),cl2,cl_stat)
    
    t_tot <- c(t_lag,t2,t_stat)
    
    
    
    df <- tibble(t_tot,cl_tot)
    
    df <- df %>% rename("Time in minutes" = t_tot) %>%
        rename("Log #cells" = cl_tot)
    
    p <- ggplot(data = df, aes(y=`Log #cells`,x=`Time in minutes`))+
        geom_line(size=1.5)+
        theme_bw()+
        theme(axis.title = element_text(size = 16, face = "bold"),
              axis.text = element_text(size = 12, face = "bold"))+
        ylab("Log cell number")+
        xlab("Time in minutes")
    
    results <- list(p,df,g)
    names(results) <- c("Growth curve", "table","G")
    return(results)
    
}



# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel(h1("Growth curve", align = "center")),
    
    sidebarLayout(
        sidebarPanel(numericInput(inputId = "gtime", label = 
                                      "Enter your calculated 
                 generation time(min)", value = NA),
                     textOutput("calcG"),
                     numericInput(inputId = "B", 
                                  label ="Calculate final cell number", 
                                  value = NA),
                     textOutput("cg"),
                     textOutput("finc")),
        mainPanel(
            #h1("Growth curve", align = "center"),
            shinycssloaders::withSpinner(plotlyOutput("gcurve"),
                                         size = 2, type = 6)
        )))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    a <- gnt3()
    clc <- c(sample(seq(60,1500,by = 50), size = 1), 
             sample(seq(100,10000,by = 100), size = 1))
    
    output$gcurve <- renderPlotly({
        Sys.sleep(0.5)
        ggplotly(a$`Growth curve`)    
        
    })
    
    output$calcG <- renderText({
        if(is.na(input$gtime)){
            "Please enter your calculated generation time \n"    
        } else {
            if(input$gtime < a$G*1.1 & input$gtime > a$G*0.9){
                print("Your generation time is correct")
            } else{
                print("Your generation time is not within 10% of the correct generation time")
            }
        }
        
    })
    
    output$cg <- renderText({
        paste("Use" , clc[1],"minutes as the time interval and", clc[2], 
              "as the initial cell number to calculate final cell number and 
          enter it to the box above.", sep=" ")
        
    })
    
    output$finc <- renderText({
        
        if(is.na(input$B)){
            "Don't forget to calculate your generation time first"
        } else{
            if(input$B < (clc[2]*2^(clc[1]/input$gtime))*1.1 
               & input$B > (clc[2]*2^(clc[1]/input$gtime))*0.9){
                "Your calculated final cell number is correct"
            } else {
                "Your calculated final cell number is incorrect"
            }
        }
        
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
