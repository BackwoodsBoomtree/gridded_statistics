library(terra)
library(stringr)

in_dir      <- "G:/ERA5/tropical_africa_ecoregions"
var_list    <- c("t2m", "vpd", "par")
out_dir    <- "G:/ERA5/tropical_africa_ecoregions_stats/"
out_name   <- "Africa_ERA5_monthly_2019-2022_"
ecos_dir    <- "G:/Africa/Tropical_Africa_Ecoregions/ecoregions"
ecos_polys  <- list.files(ecos_dir, pattern = "*.shp", full.names = TRUE, recursive = TRUE) # For naming


# if (!dir.exists(out_dir)) {
#   dir.create(out_dir, recursive = TRUE)
# }

# Remove unneeded ecoregions
ecos_polys <- ecos_polys[-c(4, 8, 9, 11)]

# Vector of ecoregion names
for (i in 1:length(ecos_polys)){
  
  eco      <- vect(ecos_polys[i])
  eco_name <- str_to_title(eco$ECO_NAME)
  eco_name <- gsub(" ", "_", eco_name)
  
  if (i == 1){
    ecos_names <- eco_name
  } else {
    ecos_names <- c(ecos_names, eco_name)
  }
}

ecos_names <- sort(ecos_names)

ncs <- list.files(in_dir, pattern = "*.nc", full.names = TRUE, recursive = TRUE)

for (i in 1:length(ncs)) {
  df        <- data.frame(matrix(ncol = 1, nrow = 48))
  
  rast_time <- rast(ncs[i], subds = var_list[1])
  dates     <- time(rast_time)
  dates     <- sub(" UTC", "", dates)
  
  df[,1]    <- dates
  colnames(df) <- "Date"

  for (j in 1:length(var_list)){
    rast_var <- rast(ncs[i], subds = var_list[j])
    m        <- global(rast_var, fun = "mean", na.rm = TRUE)$mean
    n        <- global(rast_var, fun = "notNA")$notNA
    sd       <- global(rast_var, fun = "sd", na.rm = TRUE)$sd
    sem      <- sd / sqrt(n)
    
    m_name   <- paste0(var_list[j], "_mean")
    n_name   <- paste0(var_list[j], "_n")
    sd_name  <- paste0(var_list[j], "_sd")
    sem_name <- paste0(var_list[j], "_sem")
    
    df_var   <- cbind(m, n, sd, sem)
    colnames(df_var) <- c(m_name, n_name, sd_name, sem_name)
    
    df    <- cbind(df,df_var)    
  }
  
  full_out_name <- paste0(out_dir, out_name, ecos_names[i], ".csv")
  write.csv(df, full_out_name, row.names = FALSE)
  message(paste0("Saved ", full_out_name))
  
}