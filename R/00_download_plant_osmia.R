library(tidyverse)

plants_osmia_resources <- read_delim("inputs/plants_osmia_resources.csv", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
plants_osmia_resources %>% 
  filter(`CA Native` == 'yes') %>% 
  count
