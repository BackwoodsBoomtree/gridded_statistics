library(terra)
library(stringr)
# library(plantecophys)

in_file     <- "G:/ERA5/tropical_africa/era5_tropical_africa.nc"
out_dir     <- "G:/ERA5/tropical_africa_ecoregions"
out_name    <- "/Africa_ERA5_2019-2022_"
ecos_dir    <- "G:/Africa/Tropical_Africa_Ecoregions/ecoregions"
forest_mask <- "G:/Africa/Forest_Masks/dissolved/Africa_merged_2019_2.5km_Buffer.shp"
ecos_polys  <- list.files(ecos_dir, pattern = "*.shp", full.names = TRUE, recursive = TRUE)

# Remove unneeded ecoregions
ecos_polys <- ecos_polys[-c(4, 8, 9, 11)]

# Get ERA5 data and mask by forest
era5_data   <- sds(in_file)
f_mask      <- vect(forest_mask)

# Can't mask a SpatialRasterDataset, so have to do it for each variable
d2m        <- era5_data$d2m
t2m        <- era5_data$t2m
evavt      <- era5_data$evavt
ssrd       <- era5_data$ssrd
e          <- era5_data$e
tp         <- era5_data$tp

# Convert to C, mm, and W/m2
d2m <- d2m - 273.15
t2m <- t2m - 273.15
tp  <- tp * 1000
par <- ssrd / 86400 * 0.47 # 86400 converts from J to W and 0.47 is for PAR

# Drop temps <= 0 for vpd calculation
t2m[t2m <= 0] <- NA

# Calculations from Monteith and Unsworth (2008) for temps above 0 deg C.
es       <- 0.61078 * exp((17.267 * t2m) / (237.3 + t2m)) # saturation vapor pressure or vapor pressure at air temperature
ea       <- 0.61078 * exp((17.267 * d2m) / (237.3 + d2m)) # actual vapor pressure or vapor pressure at dewpoint temperature
# rh       <- ea / es * 100
vpd <- es - ea

# Clip to forest
d2m_forest   <- mask(d2m, f_mask)
t2m_forest   <- mask(t2m, f_mask)
evavt_forest <- mask(evavt, f_mask)
par_forest   <- mask(par, f_mask)
e_forest     <- mask(e, f_mask)
tp_forest    <- mask(tp, f_mask)
vpd_forest   <- mask(vpd, f_mask)

# Fuction below can be used - result is the same as using equations
# vpd_forest   <- DewtoVPD(d2m, t2m)


for (i in 1:length(ecos_polys)){
  
  eco      <- vect(ecos_polys[i])
  eco_name <- str_to_title(eco$ECO_NAME)
  eco_name <- gsub(" ", "_", eco_name)
  
  d2m_eco   <- mask(d2m_forest, eco)
  t2m_eco   <- mask(t2m_forest, eco)
  evavt_eco <- mask(evavt_forest, eco)
  par_eco   <- mask(par_forest, eco)
  e_eco     <- mask(e_forest, eco)
  tp_eco    <- mask(tp_forest, eco)
  vpd_eco   <- mask(vpd_forest, eco)
  
  eco_data            <- sds(d2m_eco, t2m_eco, evavt_eco, par_eco, e_eco, tp_eco, vpd_eco)
  names(eco_data)     <- c(names(era5_data)[1:3], "par",
                           names(era5_data)[5:6], "vpd")
  longnames(eco_data) <- c(longnames(era5_data)[1:3], "Photosynthetically Active Radiation",
                           longnames(era5_data)[5:6], "Vapor Pressure Deficit")
  units(eco_data)     <- c("C", "C", "m of water equivalent", "W/m-2", "m of water equivalent", "mm", "kPa")
  
  out_save <- paste0(out_dir, out_name, eco_name, ".nc")
  
  writeCDF(eco_data, out_save, overwrite = TRUE, compression = 4, missval = -9999)
  message(paste0("Done! Saved ", out_save))
  
}
