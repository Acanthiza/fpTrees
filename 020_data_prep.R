

#------Setup-----

  library("tidyverse")
  library("lubridate")
  library("sf")
  source("010_data_get.R")
  

#-----Data-------

  dat <- st_read(list.files("data", full.names = TRUE, pattern = ".shp$")) %>%
    dplyr::rename(Species = SPECIES
                  , HydroGp = HYDRO_GP
                  )
  
  datTree <- dat %>%
    sf::st_set_geometry(NULL) %>%
    as_tibble %>%
    dplyr::mutate(Transect = map_chr(LOCCOMM,~substr(.,regexpr("Transect",.)[[1]],regexpr("Transect",.)[[1]]+10))
                  , Transect = as.numeric(gsub("Transect |\\.","",Transect))
                  , Year = as.numeric(format(VISITDATE,"%Y"))
                  , medYear = ((max(Year)-min(Year))/2)+min(Year)
                  , time = Year - medYear
                  , Species = as.character(Species)
                  , Species = if_else(grepl("cam",Species),"RRG",Species)
                  , Species = if_else(grepl("larg",Species),"BB",Species)
                  , Species = if_else(grepl("sten",Species),"Cooba",Species)
                  ) %>%
    dplyr::select(HydroGp,Transect,Species,Year,medYear,time,TCI,mTCI)
    
  
  