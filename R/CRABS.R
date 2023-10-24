#' CRABS: Congruent Rate Analyses in Birth-death Scenarios
#' 
#' 
#' 
#' 
#' 
#'@importFrom magrittr %>%
#'@importFrom latex2exp latex2exp
#'@importFrom tidyr gather
#'@importFrom stats approxfun rlnorm rgamma rcauchy rnorm time quantile rbeta
#'@importFrom utils txtProgressBar setTxtProgressBar head tail read.table
#'@importFrom ggplot2 ggplot aes scale_fill_gradient2 theme_classic ggtitle theme element_text scale_x_continuous 
#'@importFrom ggplot2 geom_line scale_x_reverse scale_color_manual element_blank geom_hline geom_tile xlab ylim 
#'@importFrom ggplot2 scale_fill_manual ylab labs element_rect facet_grid scale_linetype_manual scale_y_continuous geom_histogram
#'@importFrom graphics matplot par
#'@importFrom dplyr group_map bind_rows filter
#'@importFrom patchwork plot_layout
#'@importFrom colorspace sequential_hcl
#'@importFrom ape branching.times
#'@importFrom pracma fderiv
#'
#' 
#' 
#' @section References:
#' 
#' \itemize{
#'  \item Louca, S., & Pennell, M. W. (2020). Extant timetrees are consistent with a myriad of diversification histories. Nature, 580(7804), 502-505. https://doi.org/10.1038/s41586-020-2176-1
#'  \item Höhna, S., Kopperud, B. T., & Magee, A. F. (2022). CRABS: Congruent rate analyses in birth–death scenarios. Methods in Ecology and Evolution, 13, 2709– 2718. https://doi.org/10.1111/2041-210X.13997
#'  \item Kopperud, B. T., Magee, A. F., & Höhna, S. (2023). Rapidly Changing Speciation and Extinction Rates Can Be Inferred in Spite of Nonidentifiability. Proceedings of the National Academy of Sciences 120 (7): e2208851120. https://doi.org/10.1073/pnas.2208851120
#'  \item Andréoletti, J. & Morlon, H. (2023). Exploring congruent diversification histories with flexibility and parsimony. Methods in Ecology and Evolution. https://doi.org/10.1111/2041-210X.14240
#' }
#' 
"_PACKAGE"