library(terra)
# library(ncdf4)

in_file     <- "G:/ERA5/tropical_africa/era5_tropical_africa.nc"
out_dir     <- "G:ERA5/tropical_africa_ecoregions"
out_name    <- "/Africa_ERA5_2019-2022_"
ecos_dir    <- "G:/Africa/Tropical_Africa_Ecoregions/ecoregions"
forest_mask <- "G:/Africa/Forest_Masks/dissolved/Africa_merged_2019_2.5km_Buffer.shp"
ecos_polys  <- list.files(ecos_dir, pattern = "*.shp", full.names = TRUE, recursive = TRUE)

# Get ERA5 data and mask by forest
era5_data   <- sds(in_file)
f_mask      <- vect(forest_mask)

# Can't mask a SpatialRasterDataset, so have to do it for each variable
d2m        <- era5_data$d2m
d2m_forest <- mask(d2m, f_mask)

# Test ecoregion mask
d2m_egf    <- mask(d2m_forest, vect(ecos_polys[7]))

# Prepare land cover class and % filter
lc_filter <- rast("G:/MCD12C1/2020/reprocessed/percent/MCD12C1.A2020001.006.Percent_LC_03.tif")
lc_filter[lc_filter < 90] <- NA

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

list_nc <- function(in_dir, years, vi_list) {
  # Create a df of file names for each VI and year.
  # These are monthly files.
  
  ncs <- list.files(in_dir, pattern = "*.nc", full.names = TRUE, recursive = TRUE)
  
  for (i in 1:length(vi_list)) {
    
    vi_ncs <- ncs[grepl(paste0(vi_list[i], "\\."), ncs)]
    
    for (j in 1:length(years)) {
      year_ncs <- vi_ncs[grepl(years[j], vi_ncs)]
      
      if (j == 1) {
        year_out_ncs <- year_ncs
      } else {
        year_out_ncs <- c(year_out_ncs, year_ncs)
      }
    }
    
    if (i == 1) {
      df <- year_out_ncs
    } else {
      df <- cbind(df, year_out_ncs)
    }
    # year_out_ncs        <- data.frame(year_out_ncs)
  }
  colnames(df) <- vi_list
  return(df)
}
get_ts  <- function(df_f, vi_list, out_dir, roi, lc_filter) {
  
  lc_filter <- crop(lc_filter, roi, mask = TRUE, touches = TRUE)
  
  for (i in 1:length(vi_list)) {
    vi_files <- df_f[,i]
    
    for (j in 1:length(vi_files)) {
      rast_t <- rast(vi_files[j])
      rast_t <- crop(rast_t, roi, mask = TRUE, touches = TRUE)
      rast_t <- mask(rast_t, lc_filter)
      m      <- mean(values(rast_t), na.rm = TRUE)
      n      <- length(values(rast_t))
      sd     <- sd(values(rast_t), na.rm = TRUE)
      sem    <- sd / sqrt(n)
      
      row    <- cbind(m, n, sd, sem)
      
      if (j == 1) {
        df <- row
      } else {
        df <- rbind(df, row)
      }
    }
    names(df)     <- c("mean", "n", "sd", "sem")
    full_out_name <- paste0(out_dir, out_name, vi_list[i], ".csv")
    write.csv(df, full_out_name, row.names = FALSE)
    message(paste0("Saved ", full_out_name))
  }
}


df_f     <- list_nc(in_dir, years, vi_list)
get_ts(df_f, vi_list, out_dir, roi, lc_filter)
