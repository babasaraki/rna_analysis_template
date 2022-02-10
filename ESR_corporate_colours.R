#' ESR colours for corporate brand
#'
#'@example data(iris)
#'
#' ggplot(iris, aes(x = Sepal.Length, y = Sepal.Width, fill = Species))+
#'   geom_point()+
#'   scale_fill_branded_ESR()+
#'   theme_ESR()

branded_pal_ESR <- function(
  primary = "blue"
  , other = "green"
  , direction = 1
) {

  # https://www.garrickadenbuie.com/blog/custom-discrete-color-scales-for-ggplot2/
  #
  #
  ESR_branded_colors <- list(
    "blue"="#0097db"
    , "green" = "#85C659"
    , "red" = "#ec1515"
    , "yellow" = "#febf2a"
    , "purple" = "#784f96"
    , "lightgreen" ="#1abcc1"
  )

  stopifnot(primary %in% names(ESR_branded_colors))

  function(n) {
    if (n > 6) warning("Branded Colour Palette only has 6 colors.")

    if (n == 2) {
      other <- if (!other %in% names(ESR_branded_colors)) {
        other
      } else {
        ESR_branded_colors[other]
      }
      color_list <- c(other, ESR_branded_colors[primary])
    } else {
      color_list <- ESR_branded_colors[1:n]
    }

    color_list <- unname(unlist(color_list))
    if (direction >= 0) color_list else rev(color_list)
  }
}


# scale_fill_branded_ESR <- function(
#   primary = "blue"
#   , other = "lightgreen"
#   , direction = 1
# ) {
#   # https://www.garrickadenbuie.com/blog/custom-discrete-color-scales-for-ggplot2/
#   ggplot2::continuous_scale(
#     "fill", "branded"
#     , branded_pal_ESR(primary, other, direction)
#   )
# }

# scale_color_gradientn(colours = branded_pal_ESR()(2))


scale_colour_branded_ESR <- function(
  primary = "blue"
  , other = "green"
  , direction = 1
) {
  # https://www.garrickadenbuie.com/blog/custom-discrete-color-scales-for-ggplot2/
  ggplot2::discrete_scale(
    "colour", "branded"
    , branded_pal_ESR(primary, other, direction)
  )
}



