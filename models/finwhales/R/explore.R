
# plot_points <- function(occ) {
#   occ$Date_collected <- as.Date(occ$eventDate)
#   obistools::plot_map_leaflet(occ %>% filter(occ$Date_collected == as.Date('1981-07-02')))
# }
#
# plot_time <- function(occ) {
#   occ$Date_collected <- as.Date(occ$eventDate)
#   bdvis::bdcalendarheat(occ %>% filter(lubridate::year(occ$Date_collected) < 1987), "1981 - 1986")
#   bdvis::bdcalendarheat(occ %>% filter(lubridate::year(occ$Date_collected) < 1993), "1987 - 1992")
#   bdvis::bdcalendarheat(occ %>% filter(lubridate::year(occ$Date_collected) < 1999), "1993 - 1998")
#   bdvis::bdcalendarheat(occ %>% filter(lubridate::year(occ$Date_collected) < 2005), "1999 - 2004")
#   bdvis::bdcalendarheat(occ %>% filter(lubridate::year(occ$Date_collected) < 2011), "2005 - 2010")
#   bdvis::bdcalendarheat(occ %>% filter(lubridate::year(occ$Date_collected) < 2017), "2011 - 2016")
#   table(occ %>% filter(lubridate::year(occ$Date_collected) == 1981) %>% select(Date_collected))
#   obistools::plot_map(occ %>% filter(occ$Date_collected == as.Date('1981-07-02')), zoom = TRUE)
#   obistools::plot_map_leaflet(occ %>% filter(occ$Date_collected == as.Date('1981-07-02')))
# }

#' @export
fw_prepare <- function(occ, yearlimit = 1950) {
  if (!all(c("date", "month", "year") %in% names(occ))) {
    occ <- occ %>% filter(!is.na(eventDate))
    occ$date <- as.Date(occ$eventDate)
    occ$month <- lubridate::month(occ$date)
    occ$year <- lubridate::year(occ$date)
    occ <- occ %>% filter(year >= yearlimit)
    occ <- occ %>% select(decimalLongitude, decimalLatitude, date, month, year)
    colnames(occ) <- c('x', 'y', 'date', 'month', 'year')
    occ <- cbind(data.frame(occurence='presence'), occ)
  }
  occ
}

fw_plottimehist <- function(occ, field, label) {
  occ <- fw_prepare(occ)
  ggplot2::ggplot(data=occ, ggplot2::aes_string(field)) +
    ggplot2::geom_histogram(stat='count') +
    ggplot2::labs(x=label, y='Count')
}

#' @export
fw_plotyear <- function(occ) {
  fw_plottimehist(occ, 'year', 'Year')
}

#' @export
fw_plotmonth <- function(occ) {
  fw_plottimehist(occ, 'month', 'Month')
}

#' @export
fw_plotmap <- function(occ, zoom = TRUE) {
  data <- occ
  world <- borders("world", colour="gray80", fill="gray80")
  m <- ggplot() +
    world +
    geom_point(data = data, aes(x = x, y = y), size = 2, stroke = 1, alpha = 0.3, colour = "#FF368B") +
    xlab("longitude") +
    ylab("latitude")

  if (zoom) {
    xrange <- range(data$x)
    yrange <- range(data$y)
    margin <- 0.3
    dx <- margin * (xrange[2] - xrange[1])
    dy <- margin * (yrange[2] - yrange[1])
    xrange[1] <- xrange[1] - dx
    xrange[2] <- xrange[2] + dx
    yrange[1] <- yrange[1] - dy
    yrange[2] <- yrange[2] + dy
    m <- m + coord_quickmap(xlim = xrange, ylim = yrange)
  } else {
    m <- m + coord_quickmap()
  }
  return(m)
}
