# nice density plot using melted data.frame
dens_plot <- function(x, real, levels,
                     var = "variable", val = "value",
                     quantiles = c(0.025, 0.975),
                     scales = "free", ncol = 5) {


  if(!is.data.frame(real)) {
    real <- data.frame("variable" = names(real),
                       "value" = as.numeric(real))
  }

  sum_melt <- group_by(x, .data[[var]]) %>%
    summarise(
      q_low = quantile(.data[[val]], quantiles[1]),
      q_high = quantile(.data[[val]], quantiles[2]),
      mdx = median(.data[[val]]),
      quantile = if(unique(.data[[var]]) %in% real[[var]]) {
        ecdf(.data[[val]])(real[[val]][real[[var]] == unique(.data[[var]])])
      } else {
        ecdf(.data[[val]])(0)
      }
    )

  df_dens <- group_by(x, .data[[var]]) %>% do(dens = density(.[[val]]))
  df_dens <- rbind_list_base(by(df_dens, df_dens[[var]], function(x) {
    data.frame("x" = x$dens[[1]]$x,
               "y" = x$dens[[1]]$y,
               "variable" = x[[var]])
  }))

  for (v in sum_melt[[var]]) {
    df_dens <- df_dens[-which(c(df_dens[[var]] == v &
                                  (df_dens$x > sum_melt$q_high[sum_melt[[var]] == v] |
                                     df_dens$x < sum_melt$q_low[sum_melt[[var]] == v]))), ]
  }

  area_mdx <- apply(sum_melt, 1, function(x) {
    up <- which.min(abs(df_dens$x[df_dens[[var]] == x[1]] - as.numeric(x[4])))
    which(
      near(df_dens$x,df_dens$x[df_dens[[var]] == x[1]][up], tol = 0.001) &
        df_dens[[var]] == x[1]
    )
  })

  real[[var]] <- factor(real[[var]], levels = levels)
  real <- rbind(real,
                data.frame(
                  "variable"=levels[is.na(match(levels, real[[var]]))],
                  "value"=rep(0,sum(is.na(match(levels, real[[var]]))))
                ))

  g1 <- ggplot(x) +
    geom_density(aes(x = .data[[val]], y = ..density..),
                 color = "black", fill = "white") +
    geom_area(data = df_dens, mapping = aes(x = .data$x, y = .data$y),
              fill = "#002366", color = "black", inherit.aes = FALSE) +
    geom_density(aes(x = .data[[val]], y = ..density..),
                 color = "black", fill = NA) +
    geom_area(data = df_dens[unlist(area_mdx), ], mapping = aes(x = .data$x, y = .data$y),
              fill = "#002366", color = "white", inherit.aes = FALSE) +
    geom_vline(data = real, aes(xintercept = .data[[val]]),
               col = "#e62400", size = 1, linetype = "dashed") +
    facet_wrap( ~ .data[[var]], scales = scales, ncol=ncol)  + theme_bw() +
    xlab("") + ylab("Density") +
    theme(axis.text.x = element_text(),
          panel.spacing = unit(10, units = "pt"),
          axis.title.y = element_text(size = 12),
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size = 12))


}


# boostrap from the probabilities, create the densities and compare to real
coinf_plot <- function(reps = 5000, probs, levels, total, real, plot = TRUE,
                       quantiles = c(0.025, 0.975),...) {

  # sample our coinfections
  make_coinfections <- function(test_prop, ss = total) {
    return(sample(names(test_prop), size = ss, replace = TRUE, prob = test_prop))
  }

  sims <- replicate(reps, make_coinfections(probs))
  tabbed <- apply(sims, 2, function(x) {
    as.data.frame.matrix(t(table(factor(x, levels = levels))))
  })

  sim_df <- rbind_list_base(tabbed)
  sim_melt <- suppressWarnings(reshape2::melt(sim_df))

  # make the density plot
  g1 <- dens_plot(x = sim_melt, real = real,
                  levels = levels, quantiles = quantiles,
                  ...)

  if (plot) {
    print(g1)
  }
  invisible(list("plot" = g1, "vals" = sim_melt))

}


#' Reorder an x or y axis within facets
#'
#' Reorder a column before plotting with faceting, such that the values are ordered
#' within each facet. This requires two functions: \code{reorder_within} applied to
#' the column, then either \code{scale_x_reordered} or \code{scale_y_reordered} added
#' to the plot.
#' This is implemented as a bit of a hack: it appends ___ and then the facet
#' at the end of each string.
#'
#' @param x Vector to reorder.
#' @param by Vector of the same length, to use for reordering.
#' @param within Vector of the same length that will later be used for faceting
#' @param fun Function to perform within each subset to determine the resulting
#' ordering. By default, mean.
#' @param sep Separator to distinguish the two. You may want to set this manually
#' if ___ can exist within one of your labels.
#' @param ... In \code{reorder_within} arguments passed on to \code{\link{reorder}}.
#' In the scale functions, extra arguments passed on to
#' \code{\link[ggplot2]{scale_x_discrete}} or \code{\link[ggplot2]{scale_y_discrete}}.
#'
#' @source "Ordering categories within ggplot2 Facets" by Tyler Rinker:
#' \url{https://trinkerrstuff.wordpress.com/2016/12/23/ordering-categories-within-ggplot2-facets/}
#'
#' @examples
#'
#' library(tidyr)
#' library(ggplot2)
#'
#' iris_gathered <- gather(iris, metric, value, -Species)
#'
#' # reordering doesn't work within each facet (see Sepal.Width):
#' ggplot(iris_gathered, aes(reorder(Species, value), value)) +
#'   geom_boxplot() +
#'   facet_wrap(~ metric)
#'
#' # reorder_within and scale_x_reordered work.
#' # (Note that you need to set scales = "free_x" in the facet)
#' ggplot(iris_gathered, aes(reorder_within(Species, value, metric), value)) +
#'   geom_boxplot() +
#'   scale_x_reordered() +
#'   facet_wrap(~ metric, scales = "free_x")
#'
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}


#' @rdname reorder_within
scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}


#' @rdname reorder_within
scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}
