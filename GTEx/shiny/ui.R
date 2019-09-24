navbarPage("Dark Kinase Expression",
           header = singleton(includeHTML("google-analytics.js")),
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