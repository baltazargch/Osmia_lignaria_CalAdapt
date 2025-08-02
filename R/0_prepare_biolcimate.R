library(terra)
library(sf)
library(tidyverse)

# Load your custom biovar calculator using terra
source('R/udf_calculate_biovars.R')  # this should define `terra_biovars()`

#===============================
# PARAMETERS
#===============================

# General climate model setup
models <- c("INM-CM5-0", "EC-Earth3-Veg", "MIROC6", "CNRM-ESM2-1")

# Time periods and scenarios
periods <- c('2000_2020', '2021_2040', '2041_2060', '2061_2080', '2081_2100')
ssp <- c('historical', 'ssp245', 'ssp370', 'ssp585')

# Required variables
vars <- c('pr', 'tasmin', 'tasmax')

# Input/output directories
input_dir <- 'outputs/monthly_climates'
output_dir <- 'outputs/bioclim_vars'
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#===============================
# PREPARE FILE MATRIX
#===============================

# Build case table for each variable × period × scenario
table_cases <- expand_grid(vars, periods, ssp) %>% 
  mutate(
    periods = ifelse(ssp == 'historical', '2000_2014', periods),
    case = str_c(periods, '_', ssp),
    fln = str_c(input_dir, '/', vars, '_', case, '.tif')
  )

# Split into list of cases: each element corresponds to one time slice + SSP
case_list <- split(table_cases, table_cases$case)

#===============================
# CALCULATE BIOVARS
#===============================

# Loop through each case
bio_list <- map(case_list, \(case_df) {
  message("Processing: ", unique(case_df$case))
  
  # Read rasters in correct order
  fls <- case_df %>% arrange(match(vars, vars)) %>% pull(fln)
  prec <- rast(fls[1])
  tmin <- rast(fls[2])
  tmax <- rast(fls[3])
  
  # Compute biovars for each model slice
  biovars_by_model <- map(models, \(m) {
    idx <- grep(m, names(prec))
    if (length(idx) != 12) {
      warning("Model ", m, " does not have 12 months; skipping.")
      return(NULL)
    }
    terra_biovars(prec[[idx]], tmin[[idx]], tmax[[idx]])
  })
  
  names(biovars_by_model) <- models
  return(biovars_by_model)
})

#===============================
# WRITE OUTPUTS
#===============================

# Flatten and write
walk2(names(bio_list), bio_list, \(casename, model_stack) {
  walk2(names(model_stack), model_stack, \(model_name, ras) {
    if (!is.null(ras)) {
      outname <- file.path(output_dir, paste0("biovars_", casename, "_", model_name, ".tif"))
      writeRaster(ras, outname, overwrite = TRUE)
    }
  })
})
