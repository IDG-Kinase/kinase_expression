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
