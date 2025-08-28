library(sf)
library(terra)
library(tidyverse)
library(ENMeval)
library(pROC)
library(tidysdm)
library(furrr)
options(java.parameters = "-Xmx16g")

# ---- Parameters ----
set.seed(41546)
auc_min      <- 0.75                # AUC gate for "good discriminant" models
pa_sizes     <- c(2000, 4000, 7000, 9467)  # independent PA sizes
n.cores <- 14

# ---- Load plants occs and list ----
plant_list <- read_csv('outputs/records/all_species_records_and_native.csv')
natives <- plant_list %>% filter(`CA Native` == 'yes')

fls_occs <- list.files('outputs/records/plants', '.csv', full.names = T)

natives$path <- sapply(natives$species, \(x) fls_occs[grep(str_replace(x, ' ', '_'), fls_occs)], 
                       simplify = T) %>% unname()



# Feature-class set & RM grid (Merow et al. guidance)
fc_set <- c("L","LQ","LQH", "LQHT")
rm_set <- seq(0.5, 5, by = 0.5)

# Options: cloglog predictions; compute validation metrics vs partition bg
os <- list(pred.type = "cloglog",
           validation.bg = "partition",
           path = "outputs/enmeval_maxent")  # saves maxent.jar HTMLs there

# ---- Env data (historical biovars) ----
# NOTE: if your source file is already EPSG:4326, you can drop rotate/project
hist_fls <- list.files('outputs/NA_bioclim_vars/', '.tif$', full.names = T) %>% 
  stringr::str_subset('historical')

# Calculate mean current conditions
hist_rast <- rast(hist_fls) %>% split(names(.)) %>% map(mean) %>% rast()
names(hist_rast) <- paste('bio', 1:19)
hist_rast <- rotate(hist_rast)
hist_rast <- project(hist_rast, "epsg:4326")

names(hist_rast) <- c(
  "bio1"  = "Ann_Mean_Temp",
  "bio2"  = "Mean_Diurnal_Range",
  "bio3"  = "Isothermality",
  "bio4"  = "Temp_Seasonality",
  "bio5"  = "Max_Temp_Warm_Month",
  "bio6"  = "Min_Temp_Cold_Month",
  "bio7"  = "Ann_Temp_Range", 
  "bio8"  = "Mean_Temp_Wettest_Q", 
  "bio9"  = "Mean_Temp_Driest_Q",
  "bio10" = "Mean_Temp_Warm_Q",
  "bio11" = "Mean_Temp_Cold_Q",
  "bio12" = "Ann_Precip",
  "bio13" = "Precip_Wet_Month",
  "bio14" = "Precip_Dry_Month",
  "bio15" = "Precip_Season",
  "bio16" = "Precip_Wet_Q",
  "bio17" = "Precip_Dry_Q",
  "bio18" = "Precip_Warmest_Q",
  "bio19" = "Precip_Coldest_Q"
)

#start per species loop
for(r in 1:nrow(natives)){
  # r=1
  occs <- read_csv(natives$path[r], show_col_types = F)
  
  if(file.exists('outputs/models/maxent/log_models_tunning.csv')){
    spdone <- read_csv('outputs/models/maxent/log_models_tunning.csv', show_col_types = F)
    spdone <- spdone %>% filter(species == natives$species[r])
    if(nrow(spdone) >= 4) next
  }
  
  occs_sf <- st_as_sf(occs, coords=c('decimalLongitude', 'decimalLatitude'))
  st_crs(occs_sf) <- 4326
  
  dups <- terra::extract(hist_rast[[1]], occs_sf, cells=T)$cell %>% duplicated
  
  occs_filt <- occs_sf[!dups, ] %>% st_crop(hist_rast[[1]])
  occs_thin <- thin_by_dist(occs_filt, dist_min = km2m(5))
  
  # To visualize
  # plot(hist_rast[[1]])
  # plot(occs_thin$geometry, add=T)
  
  bg_aea  <- st_transform(occs_thin, 5070)
  aream <- st_union(st_buffer(bg_aea, dist = 100000)) %>% st_transform(4326) %>% 
    st_as_sf
  
  envs <- crop(hist_rast, aream, mask=T)
  
  matPreds <- as.data.frame(envs, xy = FALSE, cells = FALSE)
  varsCor  <- fuzzySim::corSelect(matPreds, var.cols = 1:ncol(matPreds))
  sel_vars <- varsCor$selected.vars
  
  stopifnot(length(sel_vars) >= 2)
  envs <-  envs[[ sel_vars ]]
  
  for(pa in pa_sizes){
    if(pa %in% spdone$npseu) next
    bg <- sample_pseudoabs(data = occs_thin, 
                           raster = envs,
                           n = pa,
                           method = "random",
                           class_label = "background",
                           return_pres = F)
    
    bg <- bg %>% filter(class == 'background') %>% st_coordinates()
    
    occs_enm <- st_coordinates(occs_thin)
    colnames(occs_enm) <- colnames(bg) <- c('longitude', 'latitude')
    
    occs_enm <- cbind(occs_enm, terra::extract(envs, occs_enm))
    bg <- cbind(bg, terra::extract(envs, bg))
    
    if(nrow(occs_enm) < 10){
      dbout <- spdone[0,]
      dbout[1,] <- NA
      dbout$species <- natives$species[1]
      dbout$noccs <- nrow(occs_enm)
      write_csv(dbout,'outputs/models/maxent/log_models_tunning.csv', append = T)
      next
    }
    n.cores.cor <- ifelse(as.numeric(st_area(aream)/1e06) > 10e06, 5, n.cores)
    
    # --- Custom evaluator: TSS on validation folds (presence vs background) ------
    e.mx <- ENMevaluate(
      occs = occs_enm[,1:2],
      envs=envs,
      bg   = bg[,1:2],
      partitions = 'randomkfold',
      tune.args   = list(fc = fc_set, rm = rm_set),
      algorithm   = "maxent.jar",
      other.settings = list(pred.type = "cloglog"),
      parallel    = TRUE, numCores = n.cores.cor,
      quiet       = FALSE
    )
    
    res <- eval.results(e.mx)%>% 
      filter(!is.na(AICc))
    if(any(res$auc.val.avg > auc_min)){
      opt_mod <-res  %>% 
        filter(auc.val.avg >= auc_min) %>% 
        filter(AICc == min(AICc, na.rm = T)) %>% 
        sample_n(1)
    } else {
      opt_mod <- res %>% 
        filter(AICc == min(AICc, na.rm = T)) %>% 
        sample_n(1)
    }
    
    if(!file.exists('outputs/models/maxent/log_models_tunning.csv')){
      dir.create('outputs/models/maxent')
      cbind(species=natives$species[r], noccs = nrow(occs_enm), 
            npseu=nrow(bg), opt_mod) %>% as_tibble %>% .[FALSE,] %>% 
        write_csv(.,'outputs/models/maxent/log_models_tunning.csv')
    }
    
    cbind(species=natives$species[r], noccs = nrow(occs_enm), 
          npseu=nrow(bg), opt_mod) %>% as_tibble %>% 
      write_csv(.,'outputs/models/maxent/log_models_tunning.csv', append = T)
    
  }
}


# ---------- helpers ----------
# features string ("LQHTP") -> maxent args
fc_args <- function(fc, rm){
  c(
    paste0("betamultiplier=", rm),
    paste0("linear=",     ifelse(grepl("L", fc), "true","false")),
    paste0("quadratic=",  ifelse(grepl("Q", fc), "true","false")),
    paste0("product=",    ifelse(grepl("P", fc), "true","false")),
    paste0("hinge=",      ifelse(grepl("H", fc), "true","false")),
    paste0("threshold=",  ifelse(grepl("T", fc), "true","false")),
    "maximumiterations=500",
    "outputformat=cloglog",     # so predictions come out as cloglog
    "responsecurves=false",
    "jackknife=false"
  )
}

read_envs <- function(paths, sel_vars){
  # path <- dbperiods$path[1]
  
  envs <- rast(paths)
  names(envs) <- paste('bio', 1:19)
  envs <- rotate(envs)
  envs <- project(envs, "epsg:4326")
  
  names(envs) <- c(
    "bio1"  = "Ann_Mean_Temp",
    "bio2"  = "Mean_Diurnal_Range",
    "bio3"  = "Isothermality",
    "bio4"  = "Temp_Seasonality",
    "bio5"  = "Max_Temp_Warm_Month",
    "bio6"  = "Min_Temp_Cold_Month",
    "bio7"  = "Ann_Temp_Range", 
    "bio8"  = "Mean_Temp_Wettest_Q", 
    "bio9"  = "Mean_Temp_Driest_Q",
    "bio10" = "Mean_Temp_Warm_Q",
    "bio11" = "Mean_Temp_Cold_Q",
    "bio12" = "Ann_Precip",
    "bio13" = "Precip_Wet_Month",
    "bio14" = "Precip_Dry_Month",
    "bio15" = "Precip_Season",
    "bio16" = "Precip_Wet_Q",
    "bio17" = "Precip_Dry_Q",
    "bio18" = "Precip_Warmest_Q",
    "bio19" = "Precip_Coldest_Q"
  )
  envs[[ sel_vars ]]
}
library(dismo)
library(glue)

# prepare future files and cases
models <- c("INM-CM5-0", "EC-Earth3-Veg", "MIROC6", "CNRM-ESM2-1")

periods <- c('2015-2044', 
             '2045-2074', 
             '2075-2100')

ssp <- c('ssp245', 'ssp370', 'ssp585')

dbperiods <- expand_grid(model = models, period = periods, ssp = ssp) %>% 
  mutate(path = str_c('outputs/NA_bioclim_vars/biovars_', period, '_', ssp, '_', model, '.tif' ), 
         case =  str_c(period, '_', ssp)) %>% 
  split(.$case)

stopifnot(all(sapply(dbperiods$path, file.exists)))

dir.create("outputs/final_maxent/models", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/final_maxent/rasters", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/final_maxent/csv",     recursive = TRUE, showWarnings = FALSE)

tune_df <- read_csv('outputs/models/maxent/log_models_tunning.csv', show_col_types = F)

best_cfg <- tune_df %>%
  group_by(species) %>%
  arrange(or.mtp.avg, desc(auc.val.avg), delta.AICc, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()
# World basemap via rnaturalearth (land polygons)
wrld_org <- rnaturalearth::ne_states(country = 'united states of america', 
                                     returnclass = "sf")

ca <- wrld_org %>% filter(name_en == 'California')

# ---------- fit, predict, save ----------
results <- best_cfg %>%
  group_by(species, npseu, fc, rm) %>%
  group_map(~{
    # .y <- best_cfg[1,]
    sp   <- .y$species
    np   <- .y$npseu
    fc1  <- .y$fc
    rm1  <- .y$rm
    
    message(glue("Fitting {sp} | npseu={np} | fc={fc1} | rm={rm1}"))
    
    occs <- read_csv(natives$path[natives$species == sp], show_col_types = F)
    
    occs_sf <- st_as_sf(occs, coords=c('decimalLongitude', 'decimalLatitude'))
    st_crs(occs_sf) <- 4326
    
    dups <- terra::extract(hist_rast[[1]], occs_sf, cells=T)$cell %>% duplicated
    
    occs_filt <- occs_sf[!dups, ] %>% st_crop(hist_rast[[1]])
    occs_thin <- thin_by_dist(occs_filt, dist_min = km2m(5))
    
    bg_aea  <- st_transform(occs_thin, 5070)
    aream <- st_union(st_buffer(bg_aea, dist = 100000)) %>% st_transform(4326) %>% 
      st_as_sf
    
    envs <- crop(hist_rast, aream, mask=T)
    
    matPreds <- as.data.frame(envs, xy = FALSE, cells = FALSE)
    varsCor  <- fuzzySim::corSelect(matPreds, var.cols = 1:ncol(matPreds))
    sel_vars <- varsCor$selected.vars
    
    stopifnot(length(sel_vars) >= 2)
    envs <-  envs[[ sel_vars ]]
    
    bg <- sample_pseudoabs(data = occs_thin, 
                           raster = envs,
                           n = np,
                           method = "random",
                           class_label = "background",
                           return_pres = F)
    
    bg <- bg %>% filter(class == 'background') %>% st_coordinates()
    
    occs_enm <- st_coordinates(occs_thin)
    colnames(occs_enm) <- colnames(bg) <- c('longitude', 'latitude')
    
    occs_enm <- cbind(occs_enm, terra::extract(envs, occs_enm))
    bg <- cbind(bg, terra::extract(envs, bg))
    
    
    # fit maxent (maxent.jar via dismo)
    mx_args <- fc_args(fc1, rm1)
    m_outdir <- glue("outputs/final_maxent/models/{sp}_fc{fc1}_rm{rm1}")
    dir.create(m_outdir, showWarnings = FALSE, recursive = TRUE)
    
    mx <- dismo::maxent(
      x    = raster::stack(envs),
      p    = as.matrix(occs_enm[, c("longitude","latitude")]),
      a    = as.matrix(bg[, c("longitude","latitude")]),
      args = mx_args,
      path = m_outdir
    )
    
    # predict (cloglog due to training args) and write
    pred_r <- predict(raster::stack(envs), mx, na.rm=T) %>% rast()          # RasterLayer (cloglog)
    
    out_pred <- glue("outputs/final_maxent/rasters/{sp}_fc{fc1}_rm{rm1}_cloglog.tif")
    terra::writeRaster(pred_r, out_pred, overwrite = TRUE)
    
    # evaluate AUC (train-on-train here; swap to a held-out set if desired)
    e <- dismo::evaluate(mx,
                         p = as.matrix(occs_enm[, c("longitude","latitude")]),
                         a = as.matrix(bg[, c("longitude","latitude")]),
                         x = envs)
    auc_train <- e@auc
    
    # MTP (minimum training presence) threshold from training presence predictions
    pres_vals <- terra::extract(pred_r, occs_enm[,1:2])
    # thr_mtp   <- min(pres_vals[,2], na.rm = TRUE)
    thr_p10 <- unname(stats::quantile(pres_vals, probs = 0.10,  # 10th percentile
                                      na.rm = TRUE, type = 7))
    
    # binary raster & write
    bin_t <- pred_r >= thr_p10
    out_bin <- glue("outputs/final_maxent/rasters/{sp}_fc{fc1}_rm{rm1}_binary_p10.tif")
    terra::writeRaster(bin_t, out_bin, overwrite = TRUE)
    
    # save model
    saveRDS(mx, file = glue("{m_outdir}/{sp}_fc{fc1}_rm{rm1}_dismo_maxent.rds"))
    
    #project to future scenarios
    plan(multisession, workers = 4)

    area_fut <- future_imap_dfc(dbperiods, \(x, n){
      
      # x <- dbperiods[[1]]
      
      binary_proj <- map(x$path, \(xpath) {
        
        newenv   <- read_envs(xpath, sel_vars) %>% 
          crop(ca, mask=T)
        
        pred_fut <- predict(raster::stack(newenv), mx, na.rm=T) %>% rast()
        # Apply the SAME threshold
        bin_fut     <- pred_fut >= thr_p10
        bin_fut
      })
      
      binary_proj <- rast(binary_proj) %>% modal(., na.rm=T) 
      
      out_bin <- glue("outputs/final_maxent/rasters/{sp}_{x$case[1]}_binary_p10.tif")
      terra::writeRaster(binary_proj, out_bin, overwrite = TRUE)
      
      area <- expanse(binary_proj, unit="km", byValue=TRUE)[2,'area'] %>% as_tibble()
      
      colnames(area) <- x$case[1]
      return(area)
    })
    
    plan(sequential)
    area_fut$species <- sp
    
    tibble(
      species = sp,
      npseu   = np,
      fc      = fc1,
      rm      = rm1,
      n_pres  = nrow(occs_enm),
      auc_train = as.numeric(auc_train),
      thr_mtp   = as.numeric(thr_mtp),
      thr_p10   = as.numeric(thr_p10),
      pred_path = out_pred,
      bin_path  = out_bin,
      model_rds = glue("{m_outdir}/{sp}_fc{fc1}_rm{rm1}_dismo_maxent.rds"), 
      current_area = expanse(crop(bin_t, ca, mask=T), unit="km", byValue=TRUE)[2,'area']
    ) %>% left_join(area_fut, by = 'species')
  }) %>% bind_rows()

readr::write_csv(results, "outputs/final_maxent/csv/final_models_metrics.csv")

results %>% 
  pivot_longer(
    current_area:last_col()
  ) %>% view
  ggplot(aes(y=value/100000, x=species))+
  geom_bar(stat='identity') +
  facet_grid(~name) +
  coord_flip()
