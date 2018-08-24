#' @export
fw_qcfilter <- function(occ, env, thindistkm = NA) {
  occ <- occ %>%
    filter(x != 0 & y != 0) %>%
    filter(x %% 0.5 != 0 | y %% 0.5 != 0)

  x <- raster::raster(env, layer=1)
  cells <- raster::cellFromXY(x, fw_xy(occ))
  occ <- occ[!duplicated(cells),]
  occ <- occ[complete.cases(raster::extract(env, fw_xy(occ))),]
  if(!is.na(thindistkm) & thindistkm > 0) {
    occ <- spThin::thin(occ, lat.col = 'y', long.col = 'x', spec.col = 'occurrence', thin.par = thindistkm)[[1]]
  }
  occ
}

#' Get random background
#' @export
fw_random_background <- function(n, mask, seed = 42) {
  bg <- cache_call(paste0('random_bg_', n, '_', seed, '_', names(mask), collapse = '_'),
                   expression({
                     set.seed(seed)
                     bg <- dismo::randomPoints(mask, n*1.1)
                     bg <- bg[complete.cases(raster::extract(mask, bg)),]
                     bg[1:n,]
                   }))
  stopifnot(nrow(bg) == n)
  bg <- cbind(data.frame(occurrence='absence'), bg)
}

#' @export
fw_write_tif <- function(x, outpath, load=TRUE) {
  tifoptions <- c("COMPRESS=DEFLATE","PREDICTOR=3", "ZLEVEL=6")
  raster::writeRaster(x, outpath, options=tifoptions, overwrite = FALSE)
  if(load) {
    raster::raster(outpath)
  } else {
    NULL
  }
}

#' @export
fw_biomod_formating_data <- function(species_xy, background_xy, species_env, background_env, eval_species_xy=NULL, eval_background_xy=NULL, eval_species_env=NULL, eval_background_env=NULL) {
  # create biomod formated data from already extracted data (useful to implement new background sampling strategies)
  resp_var <- c(rep(1, time=nrow(species_env)), rep(0, time=nrow(background_env)))
  expl_var <- rbind(species_env, background_env)
  resp_xy <- rbind(species_xy, background_xy)
  eval_resp_var <- NULL
  eval_expl_var <- NULL
  eval_resp_xy <- NULL
  if(!is.null(eval_species_env) & !is.null(eval_background_env)){
    eval_resp_var <- c(rep(1, time=nrow(eval_species_env)), rep(0, time=nrow(eval_background_env)))
    eval_expl_var <- rbind(eval_species_env, eval_background_env)
    eval_resp_xy <- rbind(eval_species_xy, eval_background_xy)
  }
  biomod2::BIOMOD_FormatingData(resp.var=resp_var
                                ,resp.xy=resp_xy
                                ,expl.var=expl_var
                                ,resp.name='Balaenoptera physalus'
                                ,eval.resp.var=eval_resp_var
                                ,eval.expl.var=eval_expl_var
                                ,eval.resp.xy=eval_resp_xy
                                ,PA.nb.rep=0) # don't extract background data (-> all data is already provided)
}

#' @export
fw_biomod_options <- function(path_to_maxent) {
  biomod2::BIOMOD_ModelingOptions(RF = list(do.classif = TRUE, ntree = 500,
                                            mtry = 'default', nodesize = 5, maxnodes = NULL),
                                  GLM = list(type = 'quadratic',
                                             interaction.level = 0, myFormula = NULL,
                                             test = 'BIC', family = 'binomial',
                                             control = glm.control(epsilon = 1e-08,maxit = 1000, trace = FALSE)),
                                  SRE = list(quant = 0.025),
                                  MAXENT.Phillips = list(path_to_maxent.jar = system.file("java", "maxent.jar", package = "finwhales"), memory_allocated = 512,
                                                         maximumiterations = 200, visible = FALSE,
                                                         linear = TRUE, quadratic = TRUE, product = FALSE,
                                                         threshold = FALSE, hinge = FALSE,
                                                         lq2lqptthreshold = 80, l2lqthreshold = 10,
                                                         hingethreshold = 15, beta_threshold = -1,
                                                         beta_categorical = -1, beta_lqp = -1,
                                                         beta_hinge = -1, defaultprevalence = 0.5))
}



#' maxsss_europe <- function(model, nocc, neurope, bg, ensemble = FALSE) {
#'   ## maximise sum of sensitivity / specificity to treshold
#'   predictions <-  biomod2::get_predictions(model, as.data.frame = TRUE)
#'
#'   is_occ_europe <- 1:NROW(predictions) %in% 1:neurope
#'   occ_fit <- predictions[is_occ_europe,,drop=FALSE]
#'
#'   ineurope <- function(lonlat) {
#'     lonlat[,1] > -34 & lonlat[,1] < 65 & lonlat[,2] > 29 & lonlat[,2] < 73
#'   }
#'   bg_fit <- predictions[(nocc+1):NROW(predictions),,drop=FALSE][ineurope(bg),,drop=FALSE]
#'
#'   if(ensemble) {
#'     occ_fit <- rowMeans(occ_fit)
#'     bg_fit <- rowMeans(bg_fit)
#'
#'     cutoff_bootstrap <- sapply(1:100, function(z){
#'       bg_fit_boot <- sample(1:NROW(bg_fit), NROW(occ_fit))
#'       stat <- biomod2::Find.Optim.Stat(Stat = "ROC",
#'                                        Fit = c(occ_fit, bg_fit[bg_fit_boot]),
#'                                        Obs = c(rep(1, NROW(occ_fit)), rep(0, NROW(occ_fit))))
#'       stat[1,"cutoff"]
#'     })
#'     median(cutoff_bootstrap)
#'   } else {
#'     tresholds <- c()
#'     for(algo_i in 1:NCOL(occ_fit)) {
#'       cutoff_bootstrap <- sapply(1:100, function(z){
#'         bg_fit_boot <- sample(1:NROW(bg_fit), NROW(occ_fit))
#'         stat <- biomod2::Find.Optim.Stat(Stat = "ROC",
#'                                          Fit = c(occ_fit[,algo_i], bg_fit[bg_fit_boot,algo_i]),
#'                                          Obs = c(rep(1, NROW(occ_fit)), rep(0, NROW(occ_fit))))
#'         stat[1,"cutoff"]
#'       })
#'       tresholds <- c(tresholds, median(cutoff_bootstrap)) ## treshold from the training data
#'     }
#'     names(tresholds) <- gsub("^.*?_RUN1_", "", colnames(occ_fit))
#'     tresholds
#'   }
#' }
#'
#' get_future_env <- function(layers, scenario="B1", year=2100) {
#'
#'   namelookup <- list(BO_sstrange="sst_range",
#'                      BO_sstmax="sst_max",
#'                      BO_sstmin="sst_min",
#'                      BO_sstmean="sst_mean",
#'                      BO_salinity="salinity")
#'   env <- NULL
#'   for(l in layers) {
#'     if(l %in% names(namelookup)) {
#'       layer <- raster(paste0("D:/a/data/BioOracle_scenarios/", scenario, "/", year, "/", namelookup[[l]], ".grd"))
#'     } else {
#'       warning("using current layer instead of future for ", l)
#'       layer <- sdmpredictors::load_layers(l, FALSE)
#'     }
#'     names(layer) <- l
#'     if(is.null(env)) {
#'       env <- layer
#'     } else {
#'       env <- stack(env, layer)
#'     }
#'   }
#'   env
#' }
#'
#' coastfilter <- function() {
#'   raster::raster("D:/a/data/coastfilter/septuple_coast_test.tif")
#' }
#'
#' #' @import raster
#' treshold_raster <- function(original, treshold, dirname=NULL) {
#'   bin <- original > treshold
#'   if(grepl("[.]tif$", original@file@name)) {
#'     outprefix <- sub("[.]tif$", paste0("_bin",round(treshold)), original@file@name)
#'   } else if(!is.null(dirname)) {
#'     outprefix <- paste0(dirname, "/", names(original), "_bin", round(treshold))
#'   }
#'
#'   bin <- write_tif(bin, paste0(outprefix, ".tif"))
#'   rc <- bin * coastfilter()
#'   write_tif(rc, paste0(outprefix, "_coast.tif"), load = FALSE)
#'   rc[values(rc) == 0] <- NA
#'   ## polygon for mapping purposes
#'   poly <- raster::rasterToPolygons(rc, dissolve = TRUE)
#'   if(!is.null(poly)) { ## e.g. when no suitable areas after climate change
#'     #poly <- raster::union(poly) ## alternatives gUnaryUnion, unionSpatialPolygons
#'     rgdal::writeOGR(obj=poly, dsn=dirname(outprefix), layer=paste0(basename(outprefix), "_coast"), driver="ESRI Shapefile")
#'   }
#'   bin
#' }
#'
#' treshold_proj <- function(proj, tresholds) {
#'   r <- NULL
#'   for(i in 1:nlayers(proj@proj@val)) {
#'     bin <- treshold_raster(raster(proj@proj@val, layer = i), tresholds[i], dirname = dirname(proj@proj@link))
#'     if(is.null(r)) {
#'       r <- bin
#'     } else {
#'       r <- stack(r, bin)
#'     }
#'   }
#'   r
#' }
#'
#' treshold_ensemble_proj <- function(ensemble_proj, treshold) {
#'   treshold_raster(ensemble_proj, treshold)
#' }
#'
#' create_ensemble_proj <- function(proj, coast) {
#'   r <- raster::mean(proj@proj@val)
#'   outpath <- sub("[.]img$", "_ensemble.tif", proj@proj@link)
#'   coastpath <- sub("[.]img$", "_ensemble_coast.tif", proj@proj@link)
#'   out <- write_tif(r, outpath)
#'   write_tif(r * coast, coastpath, load=FALSE)
#'   out
#' }


#' #' @title current_future_models
#' #' @details Build models with all data, then treshold using maxSSS with only european records
#' #'   (European range), finally predict future maps and use same treshold
#' #' @export
#' #' @import raster
#' current_future_models <- function(speciesname, layers, suffix, coast) {
#'   wd <- getwd()
#'   on.exit(setwd(wd))
#'   env <- sdmpredictors::load_layers(layers, equalarea=FALSE)
#'
#'   europe <- get_species_points(speciesname, thindistance = 30, invasive=TRUE, ineurope=TRUE, FALSE)[,1:2]
#'   europe <- europe[!duplicated(cellFromXY(env, europe)),]
#'   occurrence <- rbind(europe,
#'                       get_species_points(speciesname, thindistance = 30, invasive=TRUE, ineurope=FALSE, FALSE)[,1:2],
#'                       get_species_points(speciesname, thindistance = 30, invasive=FALSE, ineurope=FALSE, FALSE)[,1:2])
#'   occurrence <- occurrence[!duplicated(cellFromXY(env, occurrence)),]
#'   setwd("data/results")
#'
#'   if(speciesname == "Codium fragile subsp fragile") {
#'     speciesname <- "Codium fragile"
#'   }
#'   bg <- coastal_random()
#'
#'   options <- biomod_options()
#'   bg_env <- raster::extract(env, bg)
#'   #   if(NROW(occurrence) < (NROW(bg_env/2))) { ## expirement to limit number of background points
#'   #     ncenters <- min(NROW(occurrence), NROW(bg_env))
#'   #     km <- kmeans(bg_env, ncenters)
#'   #     km$cluster
#'   #     rnames <- sapply(1:ncenters, function(ki) rownames(bg_env[km$cluster==i,][which.min(dist(bg_env[km$cluster==i,], km$centers[i,])),]))
#'   #     #bg_env <- bg_env[sapply(1:ncenters)
#'   #     bg_env <- bg_env[rownames(bgenv) %in% rnames]
#'   #   }
#'   occ_env <- raster::extract(env, occurrence)
#'
#'   colnames(occurrence) <- c("x","y")
#'   colnames(bg) <- c("x","y")
#'   formatdata <- biomod_formating_data(paste0(speciesname, suffix), occurrence, bg, occ_env, bg_env)
#'   ## only training
#'   dataSplitTable <- as.matrix(data.frame(RUN1= c(rep(TRUE, nrow(occurrence)), rep(TRUE, nrow(bg)))))
#'   eval_methods <- c('KAPPA', 'TSS', 'ROC', 'FAR', 'SR', 'ACCURACY', 'BIAS', 'POD', 'CSI', 'ETS')
#'
#'   algorithms <- c("GLM", "MAXENT.Phillips", "RF", "SRE")
#'
#'   model <- biomod2::BIOMOD_Modeling(formatdata
#'                                     ,models=algorithms
#'                                     ,models.eval.meth=eval_methods
#'                                     ,models.options = options
#'                                     ,SaveObj = TRUE
#'                                     ,do.full.models = FALSE ## if true, models calibrated and evaluated with the whole dataset are done
#'                                     ,NbRunEval=1
#'                                     ,DataSplitTable=dataSplitTable
#'                                     ,VarImport=3
#'                                     ,modeling.id = gsub("BO_", "", paste0(layers, collapse= "_")))
#'
#'   tresholds <- maxsss_europe(model, nocc = NROW(occurrence), neurope = NROW(europe), bg)
#'   ensemble_treshold <- maxsss_europe(model, nocc = NROW(occurrence), neurope = NROW(europe), bg, ensemble = TRUE)
#'
#'   # clip env to europe for projection
#'   extent_eu  <- raster::extent(-40, 40, 28, 85)
#'   europe_env <- stack(crop(env, extent_eu))
#'
#'   current <- biomod2::BIOMOD_Projection(modeling.output = model, new.env = europe_env, proj.name = 'current', selected.models = 'all',
#'                                         binary.meth = NULL, compress = TRUE, build.clamping.mask = FALSE, output.format = '.img')
#'
#'
#'   current_ensemble <- create_ensemble_proj(current, coast)
#'   current_bin <- treshold_proj(current, tresholds)
#'   current_ensemble_bin <- treshold_ensemble_proj(current_ensemble, ensemble_treshold)
#'   year <- 2100
#'   for(scenario in c("A2", "A1B", "B1")) {
#'     future_env <- get_future_env(layers, scenario, year)
#'     europe_env_future <- stack(crop(future_env, extent_eu))
#'     future <- biomod2::BIOMOD_Projection(modeling.output = model, new.env = europe_env_future, proj.name = paste0("future_",scenario,"_",year), selected.models = 'all',
#'                                          binary.meth = NULL, compress = TRUE, build.clamping.mask = FALSE, output.format = '.img')
#'     future_ensemble <- create_ensemble_proj(future, coast)
#'     future_bin <- treshold_proj(future, tresholds)
#'     future_ensemble_bin <- treshold_ensemble_proj(future_ensemble, ensemble_treshold)
#'   }
#'   setwd(wd)
#' }


#' calculate_risk <- function() {
#'
#'   ## 1) create a rasterstack with all separate method and
#'   ## ensemble projections for the current and future (different scenarios)
#'   ## for sstmin, -max and -mean separately
#'   ## 2) sum
#'   ## 3) difference of sum current and sum future and convert to percentage
#'   ## for ensemble projections for the current and future (different scenarios)
#'   ## 3) mean of difference
#'   ## 4) sd of difference
#'
#' }
#'
#' #' @export
#' risk_areas_ensemble <- function() {
#'   wd <- getwd()
#'   on.exit(setwd(wd))
#'   setwd("data/results")
#'
#'   sst_layers <- c("sstmin", "sstmax", "sstmean")
#'
#'   year <- 2100
#'
#'   coast <- coastfilter()
#'   coast[values(coast) == 0] <- NA
#'
#'   for(sst in sst_layers) {
#'     current <- stack(list.files(".", paste0("proj_current_.*?[.]BO[.]",sst,"_ensemble_bin[0-9]+_coast[.]tif$"), recursive = T))
#'     current_sum <- sum(current)
#'     write_tif(current_sum, paste0("areas_risk/current_ensemble_",sst, ".tif"), load = FALSE)
#'     for(scenario in c("A2", "A1B", "B1")) {
#'       future <- stack(list.files(".", paste0("proj_future_",scenario,"_",year,"_.*?[.]BO[.]",sst,"_ensemble_bin[0-9]+_coast[.]tif$"), recursive = T))
#'       future_sum <- sum(future)
#'       write_tif(future_sum, paste0("areas_risk/future_ensemble_",sst,"_",scenario,"_",year,".tif"), load = FALSE)
#'       diff <- future_sum - current_sum
#'       write_tif(diff, paste0("areas_risk/diff_ensemble_",sst,"_",scenario,"_",year,".tif"), load = FALSE)
#'       write_tif(diff * coast, paste0("areas_risk/diff_ensemble_",sst,"_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'     }
#'   }
#'
#'   for(scenario in c("A2", "A1B", "B1")) {
#'     diffs <- stack(list.files("areas_risk", paste0("diff_ensemble_.*?_",scenario,"_",year,"[.]tif$"), full.names = TRUE))
#'     mean_diff <- write_tif(raster::calc(diffs, mean), paste0("areas_risk/mean_diff_ensemble_",scenario,"_",year,".tif"))
#'     write_tif(mean_diff*coast, paste0("areas_risk/mean_diff_ensemble_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'     sd_diff <- write_tif(raster::calc(diffs, sd), paste0("areas_risk/sd_diff_ensemble_",scenario,"_",year,".tif"))
#'     write_tif(sd_diff*coast, paste0("areas_risk/sd_diff_ensemble_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'   }
#'
#'   ## TODO same as above but for the separate methods ?
#'
#'   ## TODO calculate a turnover index (places where species dissapear, places where species go to)
#'
#'   for(sst in sst_layers) {
#'     current_files <- list.files(".", paste0("proj_current_.*?[.]BO[.]",sst,"_ensemble_bin[0-9]+_coast[.]tif$"), recursive = T)
#'     for(scenario in c("A2", "A1B", "B1")) {
#'       future_files <- list.files(".", paste0("proj_future_",scenario,"_",year,"_.*?[.]BO[.]",sst,"_ensemble_bin[0-9]+_coast[.]tif$"), recursive = T)
#'       for(i in 1:length(current_files)) {
#'         cf <- current_files[i]
#'         ff <- future_files[i]
#'         if(strsplit(cf, "/")[[1]][1] != strsplit(ff, "/")[[1]][1]) {
#'           stop("Files do not match")
#'         }
#'         sp_new <- raster(ff) - raster(cf)
#'         sp_gone <- sp_new
#'         sp_new[values(sp_new) != 1] <- 0
#'         sp_gone[values(sp_gone) != -1] <- 0
#'         if(i == 1) {
#'           sum_new <- sp_new
#'           sum_gone <- sp_gone
#'         } else {
#'           sum_new <- sum_new + sp_new
#'           sum_gone <- sum_gone + sp_gone
#'         }
#'       }
#'       write_tif(sum_new, paste0("areas_risk/sum_new_ensemble_",sst,"_",scenario,"_",year,".tif"), load = FALSE)
#'       write_tif(sum_new * coast, paste0("areas_risk/sum_new_ensemble_",sst,"_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'       write_tif(sum_gone, paste0("areas_risk/sum_gone_ensemble_",sst,"_",scenario,"_",year,".tif"), load = FALSE)
#'       write_tif(sum_gone * coast, paste0("areas_risk/sum_gone_ensemble_",sst,"_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'       turnover <- sum_new + abs(sum_gone)
#'       write_tif(turnover, paste0("areas_risk/turnover_ensemble_",sst,"_",scenario,"_",year,".tif"), load = FALSE)
#'       write_tif(turnover * coast, paste0("areas_risk/turnover_ensemble_",sst,"_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'     }
#'   }
#'
#'
#'   for(scenario in c("A2", "A1B", "B1")) {
#'     for(newgone in c("new", "gone")) {
#'       sums <- stack(list.files("areas_risk", paste0("sum_",newgone,"_ensemble_.*?_",scenario,"_",year,"[.]tif$"), full.names = TRUE))
#'       mean_sums <- write_tif(raster::calc(sums, mean), paste0("areas_risk/mean_sum_",newgone,"_ensemble_",scenario,"_",year,".tif"))
#'       write_tif(mean_sums*coast, paste0("areas_risk/mean_sum_",newgone,"_ensemble_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'       sd_sums <- write_tif(raster::calc(sums, sd), paste0("areas_risk/sd_sum_",newgone,"_ensemble_",scenario,"_",year,".tif"))
#'       write_tif(sd_sums*coast, paste0("areas_risk/sd_sum_",newgone,"_ensemble_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'     }
#'     turnovers <- stack(list.files("areas_risk", paste0("turnover_ensemble_.*?_",scenario,"_",year,"[.]tif$"), full.names = TRUE))
#'     mean_turnover <- write_tif(raster::calc(turnovers, mean), paste0("areas_risk/mean_turnover_ensemble_",scenario,"_",year,".tif"))
#'     write_tif(mean_turnover*coast, paste0("areas_risk/mean_turnover_ensemble_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'     sd_turnover <- write_tif(raster::calc(turnovers, sd), paste0("areas_risk/sd_turnover_ensemble_",scenario,"_",year,".tif"))
#'     write_tif(sd_turnover*coast, paste0("areas_risk/sd_turnover_ensemble_",scenario,"_",year,"_coast.tif"), load = FALSE)
#'   }
#'   setwd(wd)
#' }
#'
#' analyse_risk_areas <- function() {
#'   pre <- "data/results/areas_risk/"
#'   year <- 2100
#'   scenarios <- c("A2", "A1B", "B1")
#'   maps <- c("diff", "turnover")
#'
#'   corr_mean_sd <- data.frame(row.names = scenarios)
#'   for (scenario in scenarios) {
#'     for(map in maps) {
#'       m <- raster(paste0(pre, "mean_",map,"_ensemble_",scenario,"_",year,"_coast.tif"))
#'       sd <- raster(paste0(pre, "sd_",map,"_ensemble_",scenario,"_",year,"_coast.tif"))
#'       cv <- cor(abs(values(m)), values(sd), use = "complete.obs")
#'       corr_mean_sd[scenario, map] <- cv
#'     }
#'   }
#'   print(corr_mean_sd)
#'   write.csv2(corr_mean_sd, paste0(pre, "corr_mean_sd.csv"))
#'
#'   for(map in maps) {
#'     corr_m <- data.frame(row.names = scenarios)
#'     corr_s <- data.frame(row.names = scenarios)
#'     for(i in 1:(length(scenarios)-1)) {
#'       sc1 <- scenarios[i]
#'       corr_m[sc1,sc1] <- 1
#'       corr_s[sc1,sc1] <- 1
#'       m1 <- raster(paste0(pre, "mean_",map,"_ensemble_",sc1,"_",year,"_coast.tif"))
#'       s1 <- raster(paste0(pre, "sd_",map,"_ensemble_",sc1,"_",year,"_coast.tif"))
#'       for (j in (i+1):length(scenarios)) {
#'         sc2 <- scenarios[j]
#'         m2 <- raster(paste0(pre, "mean_",map,"_ensemble_",sc2,"_",year,"_coast.tif"))
#'         s2 <- raster(paste0(pre, "sd_",map,"_ensemble_",sc2,"_",year,"_coast.tif"))
#'         corr_m[sc1, sc2] <- cor(values(m1), values(m2), use = "complete.obs")
#'         corr_s[sc1, sc2] <- cor(values(s1), values(s2), use = "complete.obs")
#'         corr_m[sc2, sc1] <- corr_m[sc1, sc2]
#'         corr_s[sc2, sc1] <- corr_s[sc1, sc2]
#'         corr_m[sc2,sc2] <- 1
#'         corr_s[sc2,sc2] <- 1
#'       }
#'     }
#'     print(map)
#'     print(corr_m)
#'     print(corr_s)
#'     write.csv2(corr_m, paste0(pre, "corr_",map,"_mean_scenarios.csv"))
#'     write.csv2(corr_s, paste0(pre, "corr_",map,"_sd_scenarios.csv"))
#'   }
#'
#'
#'
#'   #   t <- stack(list.files(pre, "current_ensemble_.*?[.]tif$", full.names = TRUE))
#'   #   mt <- raster::calc(t, mean)
#'   #   mtu <- raster(paste0(pre, "mean_turnover_ensemble_B1_2100_coast.tif"))
#'   #   stu <- raster(paste0(pre, "sd_turnover_ensemble_B1_2100_coast.tif"))
#'   #
#'   #   cor(values(mt), values(mtu), use = "complete.obs")
#'   #   cor(values(mt), values(stu), use = "complete.obs")
#'   #
#'   #   m <- raster(paste0(pre, "mean_diff_ensemble_B1_2100_coast.tif"))
#'   #   sd <- raster(paste0(pre, "sd_diff_ensemble_B1_2100_coast.tif"))
#'   #   cor(values(m), values(sd), use = "complete.obs")
#' }
#'
#' # do this for all selected species
#' # write up results, create maps for all species (similar to the maps created by Eduardo
#' # (n pixel border line with striped polygon on top, use rqgis to convert raster to polygon))
#'
#'
#' # layers <- c("BO_sstrange", "BO_sstmax", "BO_parmean","BO_salinity")
#' # speciesname <- "Sargassum muticum"
#'
#' run_all_eduardo_species <- function() {
#'   coast <- coastfilter()
#'   coast[values(coast) == 0] <- NA
#'
#'   species <- c("Sargassum muticum", "Grateloupia turuturu", "Dictyota cyanoloma", "Undaria pinnatifida","Codium fragile subsp fragile")
#'   #species <- "Codium fragile subsp fragile"
#'   sstlayers <- c("BO_sstmin", "BO_sstmax", "BO_sstmean")
#'   for(speciesname in species) {
#'     for(sst in sstlayers) {
#'       layers <- c("BO_sstrange", sst, "BO_parmean","BO_salinity")
#'       current_future_models(speciesname, layers, suffix=paste0("_", sst), coast=coast)
#'     }
#'   }
#' }
#'
#' run_all_remaining_species <- function() {
#'   coast <- coastfilter()
#'   coast[values(coast) == 0] <- NA
#'
#'   sstlayers <- c("BO_sstmin", "BO_sstmax", "BO_sstmean")
#'   species <- species_norecords
#'   #   species <- c("Gracilaria vermiculophylla", "Grateloupia subpectinata",
#'   #                "Lomentaria hakodatensis", "Polysiphonia harveyi", "Polysiphonia morrowii")
#'   species <- c("Lomentaria hakodatensis", "Polysiphonia harveyi", "Polysiphonia morrowii")
#'   for(speciesname in species) {
#'     for(sst in sstlayers) {
#'       layers <- c("BO_sstrange", sst, "BO_parmean","BO_salinity")
#'       current_future_models(speciesname, layers, suffix=paste0("_", sst), coast=coast)
#'     }
#'   }
#' }
#' # sdminvasives:::run_all_eduardo_species()
#' # sdminvasives:::run_all_remaining_species()
#' # risk_areas_ensemble()
