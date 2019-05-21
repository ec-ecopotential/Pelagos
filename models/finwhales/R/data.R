#' @export
fw_obis <- function() {
  cache_file <- file.path(finwhales:::fw_checkdir('modeldata'), 'obis_all.rds')
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
fw_climate_predictors <- function() {
  options(sdmpredictors_datadir = finwhales:::fw_checkdir('modeldata/sdmpredictors'))
  layers <- c('BO2_tempmean_ss', 'BO_bathymean') # , 'BO2_carbonphytomean_ss', 'BO2_ppmean_ss', 'BO2_salinitymean_ss')

  pred <- sdmpredictors::load_layers(layers)
  chlo <- sdmpredictors::load_layers('BO2_chlomean_ss', rasterstack = FALSE)[[1]]
  chlofront <- fw_fronts(chlo@file@name)
  raster::stack(pred, chlofront)
}

#' @export
fw_future_climate <- function() {
  options(sdmpredictors_datadir = finwhales:::fw_checkdir('modeldata/sdmpredictors'))
  # layers <- c('BO2_tempmean_ss', 'BO_bathymean') # , 'BO2_carbonphytomean_ss', 'BO2_ppmean_ss', 'BO2_salinitymean_ss')
  future <- sdmpredictors::get_future_layers(c('BO2_tempmean_ss', 'BO2_chlomean_ss'), 'RCP45', 2050)

  futurepred <- sdmpredictors::load_layers(c(future$layer_code[1], 'BO_bathymean'))
  futurechlo <- sdmpredictors::load_layers(future$layer_code[2], rasterstack = FALSE)[[1]]
  futurechlofront <- raster::raster(fw_fronts(futurechlo@file@name))
  futurepred <- raster::crop(futurepred, raster::extent(futurechlofront))
  raster::stack(futurepred, futurechlofront)
}

#' @export
fw_getenv <- function(occ, ecopuserpwd, cmemsuserpwd) {
  occ$datetimes <- lubridate::as_datetime(occ$eventDate)
  occ <- mutate(occ, day = lubridate::day(datetimes), month = lubridate::month(datetimes))
  yearmonths <- occ %>% select(date_year, month) %>% distinct()
  for(ymi in 1:NROW(yearmonths)) {
    year <- yearmonths[ymi,]$date_year
    month <- yearmonths[ymi,]$month
    occfilter <- occ$date_year == year & occ$month == month
    xy <- occ[occfilter, c("decimalLongitude", "decimalLatitude")]

    env <- list(bathy = finwhales:::fw_bathymetry(),
                chlo = finwhales:::fw_fronts(finwhales:::fw_download_ecop_monthly_chlo(year, month, ecopuserpwd)),
                sst = finwhales:::fw_fronts(finwhales:::fw_download_monthly_temperature(year, month, cmemsuserpwd)))
    for(envi in seq_along(env)) {
        occ[occfilter, names(env)[envi]] <- raster::extract(env[[envi]], xy)
    }
  }
  occ[, c('bathy', 'chlo', 'sst')]
}

#' @export
fw_fronts <- function(rasterfilepath) {
  outdir <- finwhales:::fw_checkdir(paste0(dirname(rasterfilepath), '/fronts/'))
  outfile <- paste0(outdir, basename(rasterfilepath))
  if(!file.exists(outfile)) {
    x <- raster(rasterfilepath)
    fronts <- grec::detectFronts(x)
    tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=1", "ZLEVEL=6")
    raster::writeRaster(fronts, outfile, options = tifoptions, overwrite = T, NAflag = -9999)
  }
  outfile
}

fw_bathymetry <- function() {
  sdmpredictors::load_layers("BO_bathymean", equalarea = F,  rasterstack = T, datadir = "modeldata/sdmpredictors")[[1]]
}

# Starting in 2008
#' @export
fw_download_ecop_daily_sst <- function(yearmonthday, ecopuserpwd) {
  fname <- paste0(yearmonthday, "_med_sst_l4.tif")
  sstpath <- paste0("SST_L4_daily/", fname)
  outfile <- paste0(fw_checkdir('modeldata/sst/'), fname)
  fw_download_ecopotential_file(ecopuserpwd, sstpath, outfile)
  outfile
}

#' Starting in 1998
#' @export
fw_download_ecop_monthly_chlo <- function(year, month, ecopuserpwd) {
  fname <- paste0(year, stringr::str_pad(month, 2, pad = "0"), "_med_chl_l4.tif")
  chlopath <- paste0("CHL_L4_monthly/", fname)
  outfile <- paste0(fw_checkdir('modeldata/chlo/'), fname)
  fw_download_ecopotential_file(ecopuserpwd, chlopath, outfile)
  outfile
}

fw_download_ecopotential_file <- function(ecopuserpwd, ftppath, outfile) {
  rooturl <- "sftp://frontend.recas.ba.infn.it/lustre/ecopotential/incoming/PAs/LME2_Med_Pelagos_Sanctuari/EO_Biophysical/"
  if(!file.exists(outfile)) {
    print(paste0(ftppath, ' -> ', outfile))
    con <-  RCurl::getCurlHandle( ftp.use.epsv = FALSE, userpwd=ecopuserpwd)
    bin <- RCurl::getBinaryURL(paste0(rooturl, ftppath), curl = con, dirlistonly = FALSE)
    writeBin(bin, outfile)
  }
  outfile
}

#' @export
fw_download_monthly_temperature <- function(year, month, cmemsuserpwd) {
  if(year <= 2016) {
    directory <- "MEDSEA_REANALYSIS_PHYS_006_004/sv03-med-ingv-tem-rean-m/"
    outfile <- paste0('modeldata/sst/', year, stringr::str_pad(month, 2, pad = "0"), '_cmems_medsea_reanalysis.tif')
    outfile <- fw_download_copernicus_file(cmemsuserpwd, directory, outfile, year, month)
  } else if(year > 2016) {
    # TODO ftp://nrt.cmems-du.eu/Core/MEDSEA_ANALYSIS_FORECAST_PHY_006_013/sv04-med-ingv-tem-an-fc-m
    stop("TODO get data from FORECAST")
  }
  outfile
}


fw_download_copernicus_file <- function(cmemsuserpwd, directory,  outfile, year, month, day=NA) {
  rooturl <- "ftp://my.cmems-du.eu/Core/"
  if(!file.exists(outfile)) {
    # download
    # convert to compressed tif
    month <- stringr::str_pad(month, 2, pad = "0")
    ftpdir <- paste0(rooturl, directory, year, "/", month, "/")
    # ftpfile <- "20080901_mm-INGV--TEMP-MFSs4b3-MED-b20140620_re-fv05.00.nc"
    files <- RCurl::getURL(ftpdir, userpwd = cmemsuserpwd, ftp.use.epsv = FALSE, dirlistonly = TRUE)
    files <- stringr::str_split(files, '\n', simplify = TRUE)
    if(!is.na(day)) {
      ftpfile <- files[grepl(paste0('^', year, month, stringr::str_pad(day, 2, pad = "0")), files)]
      stopifnot(length(file) == 1)
    } else {
      ftpfile <- files[1]
    }
    bin <- RCurl::getBinaryURL(paste0(ftpdir, ftpfile), userpwd = cmemsuserpwd, dirlistonly = FALSE)
    writeBin(bin, paste0(outfile, ".nc"))
    b <- raster::brick(paste0(outfile, ".nc"), varname = "votemper", level = 1)
    bx <- raster::raster(b, layer = 1)
    tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=1", "ZLEVEL=6")
    raster::writeRaster(bx, outfile, options = tifoptions, overwrite = T, NAflag = -9999)
  }
  outfile
}
