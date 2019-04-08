##
# fit by ML the lambda(moi) and the p1,2,3

# function to create all moi permutations
moi_perms <- function(m, n = 3) {
  expr <- parse(text = paste0("expand.grid(",
                              paste(rep("0:m", n - 1), collapse = ","),
                              ")"))
  X <- t(as.matrix(eval(expr)))
  X <- X[, colSums(X) <= m]
  X <- t(rbind(X, m - colSums(X)))
  colnames(X) <- NULL
  return(X)

}

# density for a permutation table for given probabilities
perms_multinom_density <- function(perms, probs) {

  # first lets create all the possible combinations
  levels <- parse(text = paste0("expand.grid(",
                                paste0(rep(
                                  "c(TRUE,FALSE)", ncol(perms)
                                ), collapse = ","),
                                ")"))
  levels <- head(eval(levels), 2 ^ ncol(perms) - 1)

  # tidy the names of the probabilities
  nms <- names(probs)
  levels <-
    (apply(levels, 1, function(x) {
      paste0(nms[x], collapse = "/")
    }))
  levels <- levels[order(nchar(levels), levels)]


  # now the density of the infection composition
  size <- sum(perms[1, ])
  dens <- apply(perms, 1, dmultinom, size = size, prob = probs)

  # lastly fill these in the levels vector
  infs <- apply(perms, 1, function(x) {
    paste0(nms[x > 0], collapse = "/")
  })
  dens_inf <- vector("numeric", length = length(levels))
  names(dens_inf) <- levels

  for (i in names(dens_inf)) {
    dens_inf[i] <- sum(dens[infs == i])
  }

  return(dens_inf)
}

# create multinomial probability table for a max_moi, frequency/probailities and mean_moi
infection_probabilities <- function(max_moi, probs, lambda) {
  perms <- lapply(1:max_moi, FUN = moi_perms, n = length(probs))

  multi_dens <- lapply(perms, perms_multinom_density, probs)
  multi_dens <- do.call(rbind, multi_dens)

  w <- actuar::dztpois(1:max_moi, lambda)

  multi_dens <- apply(multi_dens, 2, weighted.mean, w)
  return(multi_dens)
}

# likelihoood of data for a given mean moi and frequency/probability of a species
ll <- function(probs, mean_moi, data, max_moi = 25) {
  probs <- infection_probabilities(max_moi = max_moi, probs, mean_moi)

  names(data) <- tolower(names(data))
  names(probs) <- tolower(names(probs))

  if (any(is.na(match(names(data), names(probs))))) {
    stop("Names in data do not match names in probs")
  }

  data_tidy <- vector(mode = "numeric", length = length(probs))
  names(data_tidy) <- names(probs)
  data_tidy[match(names(data), names(data_tidy))] <- data

  ll <- dmultinom(data_tidy, size = sum(data_tidy), prob = probs, log = TRUE)

  return(ll)
}

ll_function <- function(params, data,
                        max_moi = 25,
                        start = c("pf" = 0.7, "po" = 0.2, "pv" = 0.1, "mean_moi" = 3)) {

  # tidy the names of the params so it works in the optim
  names(params) <- names(start)

  # split the parameters as need
  probs <- params[-length(params)]
  mean_moi <- params[length(params)]
  probs <- probs / sum(probs)

  # calculate the likelihood
  ll(probs, mean_moi, data, max_moi)

}

# boostrap from the probabilities, create the densities and compare to real
coinf_plot <- function(reps = 5000, probs, levels, total, real, plot = TRUE) {

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
  sum_melt <- group_by(sim_melt, .data$variable) %>%
    summarise(
      q5 = quantile(.data$value, 0.05),
      q95 = quantile(.data$value, 0.95),
      mdx = median(.data$value),
      quantile = if(unique(.data$variable) %in% real$variable) {
        ecdf(.data$value)(real$value[real$variable == unique(.data$variable)])
         } else {
           ecdf(.data$value)(0)
         }
    )

  df_dens <- group_by(sim_melt, .data$variable) %>% do(dens = density(.$value))
  df_dens <- rbind_list_base(by(df_dens, df_dens$variable, function(x) {
    data.frame("x" = x$dens[[1]]$x,
               "y" = x$dens[[1]]$y,
               "variable" = x$variable)
  }))

  for (v in sum_melt$variable) {
    df_dens <- df_dens[-which(c(df_dens$variable == v &
                                  (df_dens$x > sum_melt$q95[sum_melt$variable == v] |
                                     df_dens$x < sum_melt$q5[sum_melt$variable == v]))), ]
  }

  area_mdx <- apply(sum_melt, 1, function(x) {
    up <- which.min(abs(df_dens$x[df_dens$variable == x[1]] - as.numeric(x[4])))
    which(
      near(df_dens$x,df_dens$x[df_dens$variable == x[1]][up], tol = 0.001) &
        df_dens$variable == x[1]
    )
  })

  real$variable <- factor(real$variable, levels = levels)
  real <- rbind(real,
                data.frame(
                  "variable"=levels[is.na(match(levels, real$variable))],
                  "value"=rep(0,sum(is.na(match(levels, real$variable))))
                ))

  g1 <- ggplot(sim_melt) +
    geom_density(aes(x = .data$value, y = ..density..),
                 color = "black", fill = "white") +
    geom_area(data = df_dens, mapping = aes(x = .data$x, y = .data$y),
              fill = "#002366", color = "black", inherit.aes = FALSE) +
    geom_density(aes(x = .data$value, y = ..density..),
                 color = "black", fill = NA) +
    geom_area(data = df_dens[unlist(area_mdx), ], mapping = aes(x = .data$x, y = .data$y),
              fill = "#002366", color = "white", inherit.aes = FALSE) +
    geom_vline(data = real, aes(xintercept = .data$value),
               col = "#e62400", size = 1, linetype = "dashed") +
    facet_wrap( ~ .data$variable, scales = "free")  + theme_bw() + xlab("") + ylab("Density") +
    theme(axis.text.x = element_text(),
      panel.spacing = unit(10, units = "pt"),
      axis.title.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 12))

  if (plot) {
  print(g1)
  }
  invisible(list("plot" = g1, "vals" = sum_melt))

}

#' Co-occurence test
#'
#' Maximum likelihood test and plotted results for coinfection dynamics
#' @param data Observed data. Either as a named vector or a 2 column
#'   \code{\link{data.frame}}. See examples for more info.
#' @param max_moi Maxuimum infction composition explored. Default = 25.
#' @param boot_iter Bootstrap iterations. Default = 10000
#' @param plot Boolean for default plotting the bootsrap results. Default = TRUE
#' @rdname cooccurence_test
#'
#' @details Estimates the maximum likely population frequency of each species
#'   (i.e. the names of the observed entities in our data) and the mean number
#'   of infections given our observed data. These estimates are used to estimate
#'   the multinomial probability distribution for all cooccurences and comparing
#'   these to our data using a boostrap method.
#'
#' @return Invisibly returns a list containing the estimated parameters
#'   and a plot of our data compared to the boostrapped estimates.
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # example of the two forms of data type accepted
#' real <- data.frame(
#' "variable"=c("pf/po","pf/pv","pf/po/pv","pf","po","po/pv","pv"),
#' "value"=c(84,179,1,5181,44,1,309)
#' )
#'
#' real <- c(
#' "pf/po" = 84, "pf/pv" = 179, "pf/po/pv" = 1,
#' "pf" = 5181, "po" = 44, "po/pv" = 1, "pv" = 309
#' )
#'
#' res <- cooccurence_test(data = real)
#'
#' }
cooccurence_test <- function(data, max_moi = 25, boot_iter = 10000, plot = TRUE) {

  # format our data
  if (is.data.frame(data)) {

    real <- data
    classes <- sapply(data, class)
    names(real) <- c("variable", "value")[c(c(1:2)[-which(classes == "numeric")],which(classes == "numeric"))]

  } else {

    real <- data.frame("variable" = names(data),
                       "value" = as.numeric(data))

  }

  real$variable <- strsplit(as.character(real$variable), "/") %>%
    lapply(sort) %>%
    lapply(paste0, collapse = "/") %>% unlist
  data <- real$value
  names(data) <- real$variable

  # useful starting parameters
  total <- sum(real$value)
  spcs <- sort(strsplit(as.character(real$variable), "/") %>% unlist %>% unique())
  start <- sapply(spcs, function(x) {
      sum(real$value[grep(x, real$variable)]) / total
    })

  # fit our best model
  fit <- optim(
    par = c(start / sum(start), 2),
    fn = ll_function,
    data = data,
    max_moi = max_moi,
    start = start,
    method = "L-BFGS-B",
    lower = c(0.0001, 0.0001, 0.0001, 0.1),
    upper = c(0.9999, 0.9999, 0.9999, max_moi),
    control = list(
      trace = TRUE,
      fnscale = -1,
      maxit = 10000
    )
  )

  # tidy the params as needed
  lambda <- tail(fit$par, 1)
  probs <- head(fit$par, length(fit$par) - 1)
  probs <- probs / sum(probs)
  names(probs) <- spcs
  fitted_multinomial <- infection_probabilities(25, probs, lambda)

  # create the bootstrapped coinf plot
  coinf <- coinf_plot(
    reps = boot_iter,
    probs = fitted_multinomial ,
    levels = names(fitted_multinomial),
    total = sum(data),
    real = real
  )

  # print the plot and return the fit and plot
  params <- list("freqs" = probs, "lambda" = lambda, "multinom" = fitted_multinomial)
  invisible(list("params" = params, "plot" = coinf))

}
