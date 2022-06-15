#' CRABS: Congruent Rate Analyses in Birth-death Scenarios
#' 
#' 
#' 
#' 
#' 
#'@importFrom magrittr %>%
#'@importFrom latex2exp latex2exp
#'@importFrom tidyr gather
#'@importFrom stats approxfun rlnorm rgamma rcauchy rnorm time
#'@importFrom utils txtProgressBar setTxtProgressBar head tail read.table
#'@importFrom ggplot2 ggplot aes scale_fill_gradient2 theme_classic ggtitle theme element_text scale_x_continuous 
#'@importFrom ggplot2 geom_line scale_x_reverse scale_color_manual element_blank geom_hline geom_tile xlab ylim 
#'@importFrom ggplot2 scale_fill_manual ylab labs element_rect facet_grid scale_linetype_manual scale_y_continuous geom_histogram
#'@importFrom dplyr group_map bind_rows filter
#'@importFrom patchwork plot_layout
#'@importFrom colorspace sequential_hcl
#'@importFrom ape branching.times
#'@importFrom pracma fderiv
#' 
#' 
#' @section References:
#' 
#' \itemize{
#'  \item Louca, S., & Pennell, M. W. (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature, 580(7804), 502-505.
#' 
#' 
#' }
#' 
"_PACKAGE"