library(shiny)
library(gganatogram)
library(tidyverse)
library(DT)
library(here)

kinase_percentiles = read_csv(here('GTEx/GTex_kinase_percentiles.csv'))
dark_kinases_percentiles = kinase_percentiles %>% 
  filter(class == "Dark")

dark_kinase_order = dark_kinases_percentiles %>%
  group_by(symbol) %>%
  summarise(num_above_90 = sum(kinase_percentile >= 90)) %>%
  arrange(desc(num_above_90))

dark_kinase_order_by_system = dark_kinases_percentiles %>%
  group_by(tissue_type) %>%
  summarise(num_above_90 = sum(kinase_percentile >= 90)) %>%
  arrange(desc(num_above_90))

GTEx_to_gganato = read_rds(here('GTEx/GTEx_tissue_to_gganato.rds'))

dataTableOptions = list(paging = FALSE, 
                        order = list(list(6, 'desc')))

ui <- navbarPage("Dark Kinase Expression",
                 tabPanel("By Kinase",
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
                                   plotOutput("anato_male_by_kinase", height="10cm",width="6cm")),
                            column(3,
                                   plotOutput("anato_female_by_kinase", height="10cm",width="6cm"))),
                          
                          # Main panel for displaying outputs ----
                          fluidRow(
                            column(12,
                                   dataTableOutput("kinase_data_summary")
                            )
                          )
                 ),
                 tabPanel("By Organ System",
                          fluidRow(
                            column(4,
                                   wellPanel(
                                     selectInput(inputId = "tissue_type",
                                                 label = "Choose an Organ System:",
                                                 choices = dark_kinase_order_by_system$tissue_type),
                                     sliderInput(inputId = "min_percentile_organ",
                                                 label = "Minimum Kinase Expression Percentile:",
                                                 min=0,max=100,value=90))),
                            column(1),
                            column(3,
                                   plotOutput("anato_male_by_organ", height="10cm",width="6cm")),
                            column(3,
                                   plotOutput("anato_female_by_organ", height="10cm",width="6cm"))
                          ),

                          # Main panel for displaying outputs ----
                          fluidRow(
                            column(12,
                                   dataTableOutput("kinase_organ_summary")
                            )
                          )
                 )

)

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  
  #############################################################################
  #Kinase Searching Functions
  #############################################################################
  this_kinase_selection <- reactive({
    dark_kinases_percentiles %>%
      filter(symbol == input$kinase, 
             kinase_percentile >= input$min_percentile)
  })
  
  output$anato_male_by_kinase <- renderPlot({
    gganatogram(data=hgMale_key %>% filter(organ %in% this_kinase_selection()$organ),
                sex = 'Male') + theme_void()
  })
  
  output$anato_female_by_kinase <- renderPlot({
    gganatogram(data=hgFemale_key %>% filter(organ %in% this_kinase_selection()$organ),
                sex = 'Female') + theme_void()
  })

  output$kinase_data_summary <- renderDataTable({
    datatable(this_kinase_selection(),
              options = dataTableOptions)
  })
  
  #############################################################################
  #Organ Searching Functions
  #############################################################################
  this_organ_selection <- reactive({
    print(input$min_percentile_organ)
    dark_kinases_percentiles %>%
      filter(tissue_type == input$tissue_type, 
             kinase_percentile >= input$min_percentile_organ)
  })
  
  output$anato_male_by_organ <- renderPlot({
    print(this_organ_selection())
    gganatogram(data=hgMale_key %>% filter(organ %in% this_organ_selection()$organ),
                sex = 'Male') + theme_void()
  })
  
  output$anato_female_by_organ <- renderPlot({
    gganatogram(data=hgFemale_key %>% filter(organ %in% this_organ_selection()$organ),
                sex = 'Female') + theme_void()
  })
  
  
  output$kinase_organ_summary <- renderDataTable({
    datatable(this_organ_selection(),
              options = dataTableOptions)
  })
  
}

shinyApp(ui, server)