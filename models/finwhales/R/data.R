#' @export
fw_obis <- function() {
  dir <- rappdirs::user_cache_dir('finwhales')
  if(!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
  cache_file <- file.path(dir, 'obis_all.rds')
  if (!file.exists(cache_file)) {
    occ <- robis::occurrence('Balaenoptera physalus')
    saveRDS(occ, cache_file)
  } else {
    occ <- readRDS(cache_file)
  }
  occ
}

#'
#' @param layers Layer codes from sdmpredictors. Default \code{c("BO2_tempmean_ss", "BO2_chlomean_ss", "BO_bathymean", "BO2_carbonphytomean_ss", "BO2_ppmean_ss", "BO2_salinitymean_ss")}
#'
#' @export
fw_climate_predictors <- function(layers = c('BO2_tempmean_ss', 'BO2_chlomean_ss', 'BO_bathymean', 'BO2_carbonphytomean_ss', 'BO2_ppmean_ss', 'BO2_salinitymean_ss')) {
  options(sdmpredictors_datadir = rappdirs::user_cache_dir('finwhales'))
  sdmpredictors::load_layers(layers)
}
