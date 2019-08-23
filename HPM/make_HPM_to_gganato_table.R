library(tibble)

GTex_to_gganato = 
  tribble(
    ~tissue_type, ~organ,
    "adult_adrenal","adrenal_gland",
    "adult_colon","colon",
    "adult_esophagus","esophagus",
    "adult_frontal_cortex","frontal_cortex",
    "adult_gallbladder","gall_bladder",
    "adult_heart","heart",
    "adult_kidney","kidney",
    "adult_liver","liver",
    "adult_lung","lung",
    "adult_ovary","ovary",
    "adult_pancreas","pancreas",
    "adult_prostate","prostate",
    "adult_rectum","rectum",
    "adult_retina", NA,
    "adult_spinal_cord","spinal_cord",
    "adult_testis","testis",
    "adult_urinary_bladder","urinary_bladder"
  ) %>%
  write_rds(here('HPM/HPM_tissue_to_gganato.rds'))