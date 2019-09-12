library(shiny)
library(gganatogram)
library(tidyverse)
library(DT)
library(here)
library(BerginskiRMisc)

###############################################################################
# Data Loading

# GTEx Data Processing ####

kinase_percentiles = read_csv(here('GTEx/GTex_kinase_percentiles.csv'))
dark_kinase_percentiles = kinase_percentiles %>% 
  filter(class == "Dark") %>% 
  #This bit of code adds a column for the positional index of the ranked kinase
  #percentile, will be used for plotting purposes
  group_by(symbol) %>% 
  arrange(desc(kinase_percentile)) %>% 
  mutate(percentile_index = 1:n()) %>% 
  ungroup() %>%
  
  #Comparible bit to add an organ percentile index
  group_by(tissue_type) %>%
  arrange(desc(kinase_percentile)) %>%
  mutate(organ_percentile_index = 1:n()) %>%
  ungroup()

dark_kinase_order = dark_kinase_percentiles %>%
  group_by(symbol) %>%
  summarise(mean_percentile = mean(kinase_percentile)) %>% 
  arrange(desc(mean_percentile))

dark_kinase_order_by_system = dark_kinase_percentiles %>%
  group_by(tissue_type) %>%
  summarise(mean_percentile = mean(kinase_percentile)) %>% 
  arrange(desc(mean_percentile))

GTEx_correlations = read_rds(here('GTEx/GTEx_expression_correlations_DK_only.rds'))

# HPM Data Processing ####

HPM_kinase_percentiles = read_rds(here('GTEx/HPM_kinase_percentiles.rds')) %>% 
  filter(class == "Dark")

HPM_dark_kinase_order = HPM_kinase_percentiles %>%
  group_by(symbol) %>%
  summarise(mean_percentile = mean(kinase_percentile)) %>% 
  arrange(desc(mean_percentile))

HPM_dark_kinase_order_by_system = HPM_kinase_percentiles %>%
  group_by(tissue_type) %>%
  summarise(mean_percentile = mean(kinase_percentile)) %>% 
  arrange(desc(mean_percentile))

# Misc Global Variables ####


dataTableOptions = list(paging = FALSE, 
                        order = list(list(6, 'desc')))

mass_spec_desc = HTML(paste0("The kinase expression percentile is calculated by from all the kinase measurements in the ",
       a(href="http://humanproteomemap.org/index.php", "Human Proteome Map"),
       ". Each kinase is then ranked within each organ system to derive kinase percentile ranking."))


###############################################################################
# Main Shiny App

ui <- navbarPage("Dark Kinase Expression",
                 header = singleton(includeScript("google-analytics.js")),
                 tabPanel("By Kinase - RNAseq",
                          # Sidebar panel for inputs ----
                          fluidRow(
                            column(3,
                                   wellPanel(
                                     selectInput(inputId = "kinase",
                                                 label = "Choose a dark kinase:",
                                                 choices = dark_kinase_order$symbol),
                                     sliderInput(inputId = "min_percentile",
                                                 label = "Minimum Kinase Expression Percentile:",
                                                 min=0,max=100,value=90)),
                                   h3("Kinase Expression Percentile"),
                                   p("The kinase expression percentile is calculated by extracting 
                                     all the kinase expression values from the GTEx Project and 
                                     then determining where each kinase ranks within each organ 
                                     system. Thus, a kinase with a kinase percentile of 90 has 
                                     an expression level in the highest 10 percent of kinases in 
                                     that organ system.")),
                                  
                            # column(1),
                            column(2,
                                   plotOutput("anato_male_by_kinase", height="10cm",width="6cm")),
                            column(2,
                                   plotOutput("anato_female_by_kinase", height="10cm",width="6cm")),
                            column(3,
                                   plotOutput("kinase_percentile_dist")),
                            column(2,
                                   h3("Most Similar Kinases"),
                                   p("Kinases whose expression profiles are most similar and associated correlation values:"),
                                   h4("Light Kinases"),
                                   htmlOutput("matchedLightKinases"),
                                   h4("Dark Kinases"),
                                   htmlOutput("matchedDarkKinases"))
                          ),
                          
                          # hr(),
                          # fluidRow(
                          #   column(4,
                          #          h3("Kinase Expression Similarity")),
                          #   column(4,
                          #          h3("Light"),
                          #          htmlOutput("matchedLightKinases")),
                          #   column(4,
                          #          h4("Dark"),
                          #          textOutput("matchedDarkKinases"))
                          # ),
                          
                          
                          hr(),
                          # Main panel for displaying outputs ----
                          fluidRow(
                            column(12,
                                   dataTableOutput("kinase_data_summary")
                            )
                          )
                 ),
                 tabPanel("By Organ System - RNAseq",
                          fluidRow(
                            column(3,
                                   wellPanel(
                                     selectInput(inputId = "tissue_type",
                                                 label = "Choose an Organ System:",
                                                 choices = dark_kinase_order_by_system$tissue_type),
                                     sliderInput(inputId = "min_percentile_organ",
                                                 label = "Minimum Kinase Expression Percentile:",
                                                 min=0,max=100,value=90)),
                                   h3("Kinase Expression Percentile"),
                                   p("The kinase expression percentile is calculated by extracting 
                                     all the kinase expression values from the GTEx Project and 
                                     then determining where each kinase ranks within each organ 
                                     system. Thus, a kinase with a kinase percentile of 90 has 
                                     an expression level in the highest 10 percent of kinases in 
                                     that organ system.")),
                            column(1),
                            column(2,
                                   plotOutput("anato_male_by_organ", height="10cm",width="6cm")),
                            column(2,
                                   plotOutput("anato_female_by_organ", height="10cm",width="6cm")),
                            column(4,
                                   plotOutput("organ_percentile_dist"))
                            
                          ),
                          hr(),
                          # Main panel for displaying outputs ----
                          fluidRow(
                            column(12,
                                   dataTableOutput("kinase_organ_summary")
                            )
                          )
                 ),
                 tabPanel("By Kinase - Mass Spec",
                          # Sidebar panel for inputs ----
                          fluidRow(
                            column(3,
                                   wellPanel(
                                     selectInput(inputId = "HPM_kinase",
                                                 label = "Choose a dark kinase:",
                                                 choices = HPM_dark_kinase_order$symbol),
                                     sliderInput(inputId = "HPM_min_percentile_kinase",
                                                 label = "Minimum Kinase Expression Percentile:",
                                                 min=0,max=100,value=90)),
                                   h3("Kinase Mass Spec Level Percentile"),
                                   p(mass_spec_desc)),
                            column(1),
                            column(2,
                                   plotOutput("HPM_anato_male_by_kinase", height="10cm",width="6cm")),
                            column(2,
                                   plotOutput("HPM_anato_female_by_kinase", height="10cm",width="6cm")),
                            column(4,
                                   plotOutput("HPM_kinase_percentile_dist"))
                          ),
                          hr(),
                          # Main panel for displaying outputs ----
                          fluidRow(
                            column(12,
                                   dataTableOutput("HPM_kinase_data_summary")
                            )
                          )
                 ),
                 tabPanel("By Organ System - Mass Spec",
                          fluidRow(
                            column(3,
                                   wellPanel(
                                     selectInput(inputId = "HPM_tissue_type",
                                                 label = "Choose an Organ System:",
                                                 choices = HPM_dark_kinase_order_by_system$tissue_type),
                                     sliderInput(inputId = "HPM_min_percentile_organ",
                                                 label = "Minimum Kinase Expression Percentile:",
                                                 min=0,max=100,value=90)),
                                   h3("Kinase Mass Spec Level Percentile"),
                                   p(mass_spec_desc)),
                            column(1),
                            column(2,
                                   plotOutput("HPM_anato_male_by_organ", height="10cm",width="6cm")),
                            column(2,
                                   plotOutput("HPM_anato_female_by_organ", height="10cm",width="6cm")),
                            column(4,
                                   plotOutput("HPM_organ_percentile_dist"))
                            
                          ),
                          hr(),
                          # Main panel for displaying outputs ----
                          fluidRow(
                            column(12,
                                   dataTableOutput("HPM_kinase_organ_summary")
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
    dark_kinase_percentiles %>%
      filter(symbol == input$kinase, 
             kinase_percentile >= input$min_percentile) %>% 
      arrange(desc(kinase_percentile), desc(TPM))
  })
  
  output$anato_male_by_kinase <- renderPlot({
    gganatogram(data=hgMale_key %>% filter(organ %in% this_kinase_selection()$organ),
                sex = 'Male') + theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0, "cm"))
  })
  
  output$anato_female_by_kinase <- renderPlot({
    gganatogram(data=hgFemale_key %>% filter(organ %in% this_kinase_selection()$organ),
                sex = 'Female') + theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0, "cm"))
  })
  
  output$kinase_percentile_dist <- renderPlot({
    full_kinase_data = dark_kinase_percentiles %>%
      filter(symbol == input$kinase)
    
    ggplot(dark_kinase_percentiles, aes(x=percentile_index,y=kinase_percentile, group=symbol)) + 
      geom_line(alpha=0.05) + 
      ylim(c(0,100)) + 
      geom_line(data=full_kinase_data, aes(x=percentile_index,y=kinase_percentile)) +
      geom_hline(yintercept = input$min_percentile, color='red', alpha=0.5) + 
      geom_text(x=0, y=input$min_percentile, 
                label = "Min Kinase Percentile", 
                color='red', hjust=0, vjust=-0.25, alpha=0.5) +
      labs(x="Organ Index", y="Kinase Percentile") +
      theme(text = element_text(size=20),
            plot.margin = margin(0, 1, 0, 0, "cm")) +
      theme_berginski()
  })
  
  output$kinase_data_summary <- renderDataTable({
    datatable(this_kinase_selection() %>%
                select(-percentile_index,-organ_percentile_index),
              options = dataTableOptions)
  })
  
  output$matchedLightKinases <- renderPrint({
    matched_light = GTEx_correlations %>%
      filter(kinase_1 == input$kinase, kinase_2_class == "Light") %>%
      top_n(5,kinase_percentile_cor) %>%
      arrange(desc(kinase_percentile_cor))
    
    tags$ul(
      tags$li(matched_light$kinase_2[1], "-", matched_light$kinase_percentile_cor[1]),
      tags$li(matched_light$kinase_2[2], "-", matched_light$kinase_percentile_cor[2]),
      tags$li(matched_light$kinase_2[3], "-", matched_light$kinase_percentile_cor[3]),
      tags$li(matched_light$kinase_2[4], "-", matched_light$kinase_percentile_cor[4]),
      tags$li(matched_light$kinase_2[5], "-", matched_light$kinase_percentile_cor[5])
    )
  })
  
  output$matchedDarkKinases <- renderPrint({
    matched_light = GTEx_correlations %>%
      filter(kinase_1 == input$kinase, kinase_2_class == "Dark") %>%
      top_n(5,kinase_percentile_cor) %>%
      arrange(desc(kinase_percentile_cor))
    
    tags$ul(
      tags$li(matched_light$kinase_2[1], "-", matched_light$kinase_percentile_cor[1]),
      tags$li(matched_light$kinase_2[2], "-", matched_light$kinase_percentile_cor[2]),
      tags$li(matched_light$kinase_2[3], "-", matched_light$kinase_percentile_cor[3]),
      tags$li(matched_light$kinase_2[4], "-", matched_light$kinase_percentile_cor[4]),
      tags$li(matched_light$kinase_2[5], "-", matched_light$kinase_percentile_cor[5])
    )
  })
  
  
  #############################################################################
  #Organ Searching Functions
  #############################################################################
  this_organ_selection <- reactive({
    dark_kinase_percentiles %>%
      filter(tissue_type == input$tissue_type, 
             kinase_percentile >= input$min_percentile_organ) %>% 
      arrange(desc(kinase_percentile), desc(TPM))
  })
  
  output$anato_male_by_organ <- renderPlot({
    gganatogram(data=hgMale_key %>% filter(organ %in% this_organ_selection()$organ),
                sex = 'Male') + theme_void()
  })
  
  output$anato_female_by_organ <- renderPlot({
    gganatogram(data=hgFemale_key %>% filter(organ %in% this_organ_selection()$organ),
                sex = 'Female') + theme_void()
  })
  
  output$organ_percentile_dist <- renderPlot({
    full_organ_data = dark_kinase_percentiles %>%
      filter(tissue_type == input$tissue_type)
    
    ggplot(dark_kinase_percentiles, aes(x=organ_percentile_index,y=kinase_percentile,group=tissue_type)) + 
      geom_line(alpha=0.05) + 
      ylim(c(0,100)) + 
      geom_line(data=full_organ_data, aes(x=organ_percentile_index,y=kinase_percentile)) +
      geom_hline(yintercept = input$min_percentile_organ, color='red', alpha=0.5) + 
      geom_text(x=0, y=input$min_percentile_organ, 
                label = "Min Kinase Percentile", 
                color='red', hjust=0, vjust=-0.25, alpha=0.5) +
      labs(x="Dark Kinase Index", y="Kinase Percentile") +
      theme(text = element_text(size=20)) +
      theme_berginski()
  })
  
  output$kinase_organ_summary <- renderDataTable({
    datatable(this_organ_selection() %>%
                select(-percentile_index,-organ_percentile_index),
              options = dataTableOptions)
  })
  
  
  #############################################################################
  #HPM - Kinase Searching Functions
  #############################################################################
  this_HPM_kinase_selection <- reactive({
    HPM_kinase_percentiles %>%
      filter(symbol == input$HPM_kinase, 
             kinase_percentile >= input$HPM_min_percentile_kinase) %>% 
      arrange(desc(kinase_percentile), desc(total_protein_level))
  })
  
  output$HPM_anato_male_by_kinase <- renderPlot({
    gganatogram(data=hgMale_key %>% filter(organ %in% this_HPM_kinase_selection()$organ),
                sex = 'Male') + theme_void()
  })
  
  output$HPM_anato_female_by_kinase <- renderPlot({
    gganatogram(data=hgFemale_key %>% filter(organ %in% this_HPM_kinase_selection()$organ),
                sex = 'Female') + theme_void()
  })
  
  output$HPM_kinase_percentile_dist <- renderPlot({
    full_kinase_data = HPM_kinase_percentiles %>%
      filter(symbol == input$HPM_kinase)
    
    ggplot(HPM_kinase_percentiles, aes(x=percentile_index,y=kinase_percentile, group=symbol)) + 
      geom_line(alpha=0.05) + 
      ylim(c(0,100)) + 
      geom_line(data=full_kinase_data, aes(x=percentile_index,y=kinase_percentile)) +
      geom_hline(yintercept = input$HPM_min_percentile_kinase, color='red', alpha=0.5) + 
      geom_text(x=1, y=input$min_percentile, 
                label = "Min Kinase Percentile", 
                color='red', hjust=0, vjust=-0.25, alpha=0.5) +
      labs(x="Organ Index", y="Kinase Percentile") +
      theme(text = element_text(size=20)) +
      theme_berginski()
  })
  
  output$HPM_kinase_data_summary <- renderDataTable({
    datatable(this_HPM_kinase_selection() %>%
                select(-percentile_index,-organ_percentile_index),
              options = dataTableOptions)
  })
  
  
  
  
  #############################################################################
  #HPM - Organ Searching Functions
  #############################################################################
  this_HPM_organ_selection <- reactive({
    HPM_kinase_percentiles %>%
      filter(tissue_type == input$HPM_tissue_type, 
             kinase_percentile >= input$HPM_min_percentile_organ) %>% 
      arrange(desc(kinase_percentile), desc(total_protein_level))
  })
  
  output$HPM_anato_male_by_organ <- renderPlot({
    gganatogram(data=hgMale_key %>% filter(organ %in% this_HPM_organ_selection()$organ),
                sex = 'Male') + theme_void()
  })
  
  output$HPM_anato_female_by_organ <- renderPlot({
    gganatogram(data=hgFemale_key %>% filter(organ %in% this_HPM_organ_selection()$organ),
                sex = 'Female') + theme_void()
  })
  
  output$HPM_organ_percentile_dist <- renderPlot({
    full_organ_data = HPM_kinase_percentiles %>%
      filter(tissue_type == input$HPM_tissue_type)
    
    ggplot(HPM_kinase_percentiles, aes(x=organ_percentile_index,y=kinase_percentile,group=tissue_type)) + 
      geom_line(alpha=0.05) + 
      ylim(c(0,100)) + 
      geom_line(data=full_organ_data, aes(x=organ_percentile_index,y=kinase_percentile)) +
      geom_hline(yintercept = input$HPM_min_percentile_organ, color='red', alpha=0.5) + 
      geom_text(x=1, y=input$min_percentile_organ, 
                label = "Min Kinase Percentile", 
                color='red', hjust=0, vjust=-0.25, alpha=0.5) +
      labs(x="Dark Kinase Index", y="Kinase Percentile") +
      theme(text = element_text(size=20)) +
      theme_berginski()
  })
  
  output$HPM_kinase_organ_summary <- renderDataTable({
    datatable(this_HPM_organ_selection() %>%
                select(-percentile_index,-organ_percentile_index),
              options = dataTableOptions)
  })
  
  
}

shinyApp(ui, server)