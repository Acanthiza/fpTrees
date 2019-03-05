
  library("tibble")
  library("magrittr")
  library("purrr")

  fs <- list.files("R:/DFW_CBD/Projects/Projects_Ecology/RRP_2018/DataForAssessments/Katarapko", full.names = TRUE)
  
  fs <- grep("lock",fs,value = TRUE, invert = TRUE)
  
  lapply(fs,function(x) file.copy(x,"data",overwrite = TRUE))
  