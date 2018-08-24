#' @export
fw_xy <- function(df) {
  if(all(c('x', 'y') %in% colnames(df))) {
    return(df[,c('x', 'y')])
  }
  stop('x and or y column not found in df')
}

#' @export
fw_true <- function(x) {
  !is.na(x) & x
}

#' @export
fw_false <- function(x) {
  !is.na(x) & !x
}

list_cache <- function() {
  list.files(rappdirs::user_cache_dir("finwhales"), "call_", full.names = TRUE)
}

clear_cache <- function(age=36) {
  cachefiles <- list_cache()
  rmfiles <- cachefiles[difftime(Sys.time(), file.info(cachefiles)[,"mtime"], units = "hours") > age]
  unlink(rmfiles)
}

cache_call <- function(key, expr, env = NULL) {
  stopifnot(is.expression(expr))
  if(is.null(env)) {
    env = parent.frame()
  }
  cache_dir <- rappdirs::user_cache_dir("finwhales")
  cachefile <- file.path(cache_dir, paste0("call_", digest::digest(list(key=key, expr=expr)), ".rds"))
  if(file.exists(cachefile) && difftime(Sys.time(), file.info(cachefile)[,"mtime"], units = "hours") < 36) {
    return(readRDS(cachefile))
  } else {
    result <- eval(expr, envir = NULL, enclos = env)
    if(!dir.exists(cache_dir)) {
      dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    }
    saveRDS(result, cachefile)
    return(result)
  }
}
