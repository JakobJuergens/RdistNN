# Shared structural inspection for plotting contracts.
.plot_layer_geoms <- function(p) {
  vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
}
