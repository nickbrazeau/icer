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
