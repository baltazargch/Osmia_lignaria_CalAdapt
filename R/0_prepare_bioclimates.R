library(terra)
library(sf)
library(tidyverse)
library(ncdf4)
library(foreach)
library(doParallel)

# Function to aggregate each model raster to 12 monthly means
collapse_to_monthly_climatology <- function(r_model) {
  # Extract month names from layer names
  month_labels <- str_extract(names(r_model), "(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)")
  
  monthly_means <- map(month.abb, function(mon) {
    idx <- which(month_labels == mon)
    if (length(idx) > 0) {
      mean(r_model[[idx]], na.rm = TRUE)
    } else {
      NA
    }
  })
  
  r_monthly <- rast(monthly_means)
  names(r_monthly) <- month.abb
  return(r_monthly)
}

# Target models and parameters
modelsIN = c("INM-CM5-0", "EC-Earth3-Veg", "MIROC6", "CNRM-ESM2-1")

# get files names and prepare output directories
clim.files <- list.files('inputs/climates', '.nc$', full.names=T) 

# # get baseline time and models 
# hist <- clim.files[ grep('historical', clim.files) ]

dir.create('outputs/monthly_climates', recursive = TRUE, showWarnings = FALSE)

# set parallel tunning
n_cores <- 6
cl <- makeCluster(n_cores)
registerDoParallel(cl)

foreach(i = seq_along(clim.files), 
        .packages = c('terra', 'sf', 'tidyverse', 'ncdf4')) %dopar% {
          # for(i in seq_along(clim.files)) {
          # for(i in 1:2) {
          # i = 2
          # load rasters
          r <- rast(clim.files[i])
          
          # Determine variable name from filename
          out.name.var <- case_when(
            str_detect(clim.files[i], "pr_") ~ "pr",
            str_detect(clim.files[i], "tasmax_") ~ "tasmax",
            TRUE ~ "tasmin"
          )
          
          # Get actual year range from raster time dimension
          year.range <- time(r) %>%
            as_date(origin = "2000-01-01") %>%  # If time is in numeric days since 2000-01-01
            year() %>%
            range()
          
          # Combine var + year range (e.g., "tasmax_2000_2020")
          out.name.var <- paste0(out.name.var, "_", str_c(year.range, collapse = "_"), 
                                 '_', str_extract(clim.files[i], "(ssp\\d+|historical)") )
          
          if(file.exists(paste0('outputs/monthly_climates/', out.name.var, '.tif'))) {
            cat(paste0(out.name.var, ' already saved\n'))
            return(NULL)
          }
          
          num.years <- diff(year.range)+ 1
          
          nc <- nc_open(clim.files[i])
          
          model.names <- ncvar_get( nc, 'simulation' )
          
          stopifnot(nlyr(r) == length(model.names) * num.years * 12)
          
          names.vec <- rep(model.names, each = num.years * 12)
          
          times.vec <- time(r) %>%
            month(label = TRUE, abbr = TRUE) %>%
            as.character()
          
          names(r) <- str_c(names.vec, times.vec, sep = '_') 
          
          # Split model and month from layer names
          layer_names <- names(r)
          month_names <- str_extract(layer_names, "(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)")
          model_names <- str_remove(layer_names, "_(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$")
          
          # Check structure
          # table(model_names, month_names)
          
          # Create list of rasters per model
          model_list <- modelsIN
          r.list <- setNames(vector("list", length(model_list)), model_list)
          
          cat(' Splitting by model name\n')
          for (m in model_list) {
            r.list[[m]] <- r[[grep(m, names(r))]]
          }
          
          cat('Calculating monthly means\n')
          # Apply to all models
          
          # r.monthly.by.model <- map(r.list, collapse_to_monthly_climatology)
          r.monthly.by.model <- map(r.list, collapse_to_monthly_climatology)
          
          
          cat('Combining means and models\n')
          # Combine all into a single stack
          r_all <- do.call(c, r.monthly.by.model)
          
          
          cat('Writing means and models\n')
          writeRaster(rast(r_all), paste0('outputs/monthly_climates/', out.name.var, '.tif'))
        }
stopCluster(cl)
