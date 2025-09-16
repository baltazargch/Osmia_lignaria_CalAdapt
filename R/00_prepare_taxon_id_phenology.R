library(tidyverse)

plants <- read_csv('outputs/records/all_species_records_and_native.csv') %>% 
  rename_with(janitor::make_clean_names) %>% 
  filter(ca_native == 'yes')



cch2_taxonid <- read_csv('inputs/plants_phenology/Native Plants of California_1758044999.csv')

cch2_in <- cch2_taxonid %>% filter(ScientificName  %in% plants$species)

plants$species[!plants$species  %in% cch2_in$ScientificName]
