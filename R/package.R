#' \pkg{icer} Infection compistion estimation
#'
#'
#' @docType package
#' @name icer
#'
#' @importFrom stats density dmultinom ecdf median optim quantile weighted.mean
#' @importFrom stats dnbinom dpois
#' @importFrom utils head tail
#' @importFrom dplyr group_by near summarise do
#' @importFrom ggplot2 geom_density geom_vline ggplot theme theme_bw unit
#' @importFrom ggplot2 aes element_text facet_wrap geom_area xlab ylab
#' @importFrom rlang .data
#'
"_PACKAGE"

if(getRversion() >= "2.15.1")  utils::globalVariables(c(".","..density.."))
