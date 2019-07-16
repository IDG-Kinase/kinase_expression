library(shiny)
library(gganatogram)
library(tidyverse)
library(DT)
library(here)

kinase_percentiles = read_csv(here('GTEx/GTex_kinase_percentiles.csv'))
dark_kinases = kinase_percentiles %>% 
  filter(class == "Dark")

dark_kinase_order = dark_kinases %>%
  group_by(symbol) %>%
  summarise(num_above_90 = sum(kinase_percentile >= 90)) %>%
  arrange(desc(num_above_90))

ui <- fluidPage(
  # App title ----
  headerPanel("Dark Kinase Expression"),
  
  # Sidebar panel for inputs ----
  fluidRow(
    column(4,
           wellPanel(
             selectInput(inputId = "kinase",
                           label = "Choose a dark kinase:",
                           choices = dark_kinase_order$symbol),
               sliderInput(inputId = "min_percentile",
                           label = "Minimum Kinase Expression Percentile:",
                           min=0,max=100,value=90))),
    column(1),
    column(3,
           plotOutput("kinase_anato_male", height="10cm",width="6cm")),
    column(3,
           plotOutput("kinase_anato_female", height="10cm",width="6cm"))),
  
  # Main panel for displaying outputs ----
  fluidRow(
    column(12,
    dataTableOutput("kinase_data_summary")
    )
  )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  this_kinase_selection <- reactive({dark_kinases %>%
    filter(symbol == input$kinase, kinase_percentile >= input$min_percentile)
  })
  
  output$kinase_anato_male <- renderPlot({
    gganatogram(data=hgMale_key %>% filter(organ %in% this_kinase_selection()$organ),
                sex = 'Male') + theme_void()
  })
  
  output$kinase_anato_female <- renderPlot({
    gganatogram(data=hgFemale_key %>% filter(organ %in% this_kinase_selection()$organ),
                sex = 'Female') + theme_void()
  })
  
  output$kinase_data_summary <- renderDataTable({
    this_kinase_selection()
  })
  
}

shinyApp(ui, server)