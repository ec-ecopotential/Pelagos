plot_densities <- function(occ_env_train, occ_env_test, bg_env_train, bg_env_test) {

  p <- ggplot(dat, aes(x = layer, fill = Legend, colour = Legend)) +
    geom_density(alpha = 0.5) +
    labs(list(x=layers_lookup[[layer]], fill="", colour = "")) +
    scale_fill_manual(values=cbPalette[c(4,1,2)]) +
    scale_colour_manual(values=cbPalette[c(4,1,2)]) +
    facet_grid(. ~ split)
}
