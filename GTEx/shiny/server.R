function(input, output, session) {
  
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (!is.null(query[['input_kinase']]) && 
        query[['input_kinase']] %in% dark_kinase_percentiles$symbol) {
      updateSelectInput(session, "kinase", selected = query[['input_kinase']])
    }
  })
  
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