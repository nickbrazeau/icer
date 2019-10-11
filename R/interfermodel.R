# TODO need to make interfere work within new framework

#' Interference mutlinomial distribution
#'
#' @param x Vector of length equal to `data` of integers in `0:sum(data)`.
#' @param prob Numeric non-negative vector of length `length(data)`, specifying
#'   the probability for classes; is internally normalized to sum 1.
#'   Infinite and missing values are not allowed.
#' @param k Interference parameters. See details.
#' @param interf Interference relationships. Output of \code{interf_pair_pos}
#'
#' @details Interference parameters describe how the probability density of
#'   `dmultionom(x, prob=prob)` is altered depending on which is the first draw
#'   from the mutlinomial distribution.
#'
#'   Two types of interference parameter can be specified. For example:
#'
#'   Modulators, which should be passed as a named vector, with names "k<n>":
#'
#'   `k = c("k0"=0.5)` - The probability of choosing the same event as the first
#'    drawn event will be mutliplied by 0.5
#'   `k = c("k1"=0.5,"k2"=2)` - The probability of choosing the first element in
#'    x again if it was drawn first will be multiplied by 0.5, whereas if the
#'    second element was drawn first then the probability of picking it again is
#'    doubled.
#'    drawn event will be mutliplied by 0.5
#'
#'   Interferences, passed as a named vector, with names "k_<n>_<m>":
#'
#'   `k = c("k_1_2" = 0.5, "k_1_3" = 2, "k_2_3")` - The probability of choosing
#'    event 1 if event 2 drawn first is half, and vice versa. The probability
#'    of drawing event 1 if event 3 was drawn first is doubled, and vice versa.
#'    The probability of drawing event 2 if event 3 was drawn first is doubled,
#'    and vice versa.
#'    drawn event will be mutliplied by 0.5
#'    Interferences are thus symmetrical between events.
#'
#' @return Returns probability density vector of length equal to `length(prob)`
interf_dmultinom <- function(x, prob, k, interf){

  # which positions in x could have been drawn first
  pos <- which(x > 0)

  # create our density to be added to for each possible first draw
  dens <- 0

  # if interference
  if(!is.logical(interf)){

    # loop through each positon in x that could be decreased
    for(i in pos){
      x_new <- x
      x_new[i] <- x_new[i]-1
      prob_new <- prob
      prob_new[interf$p[[i]]] <- prob_new[interf$p[[i]]]*k[interf$k[[i]]]
      dens <- dens + dmultinom(x_new, prob=prob_new) * prob[i]
    }

    # otherwise we are doing modulation
  } else {

    for(i in pos){
      x_new <- x
      x_new[i] <- x_new[i]-1
      prob_new <- prob
      prob_new[i] <- prob_new[i]*k[i]
      dens <- dens + dmultinom(x_new,prob=prob_new) * prob[i]
    }

  }
  return(dens)
}

#' @noRd
interf_pair_pos <- function(k, n) {

  if(all(grepl("_", names(k)))){
    interf <- list()
    interf$k <- lapply(seq_len(n), grep, names(k))
    interf$p <- lapply(seq_len(n), function(x){
      p_ch <- grep(x, names(k))
      v <- unlist(strsplit(gsub("k_([0-9]+)$", "\\1", names(k[p_ch])),""))
      v <- v[v!=as.numeric(x)]
      return(as.numeric(v))
    })
  } else {
    interf <- FALSE
  }

  return(interf)
}

#' Interference model of strain acquisition
#'
#' @param params Vector of parameters including the frequency of each species
#'   and the mean moi mu.
#' @param data Real data to compare model predictions to.
#' @param pci_list Permutation, composition and infections list.
#' @param max_moi Maxuimum infction composition explored. Default = 25.
#'
#' @details For a given frequency of each strain and mean moi, returns the
#'   related multionomial distribution describing each infection type. `params`
#'   vector can contain extra parameters `k...` that describe the interference.
#'   See \code{interf_dmultinom} for how to format these parameters.
#'
#' @return Returns a named multinomial probability distribution.
interference <- function(params, data, pci_list,  max_moi = 25) {

  # handle frequency params
  # ------------------------------------
  n_sp <- ncol(pci_list$perms[[1]])
  freqs <- params[seq_len(n_sp)]
  freqs <- freqs / sum(freqs)

  # grab parameters related to interference
  k <- params[grep("k_|k\\d", names(params))]
  if(length(k) == 0){
    k <- 1
  }

  # adjust k if it is only one parameter long
  if(length(k) == 1){
    k <- rep(k, n_sp)
  }

  # what is the interference relationship
  interf <- interf_pair_pos(k, n_sp)

  # Density of each infection composition using a multinomial distribution
  densities <- lapply(pci_list$perms, function(x){
    dens <- apply(x, 1, interf_dmultinom, prob = freqs, k = k, interf = interf)
  })

  # turn these into a a matrix of multinomial densities
  multi_dens <- comp_dens_to_multi_dens(pci_list = pci_list, densities = densities)

  # calculate the overall multinomial distributon
  probs <- infection_probabilities(multi_dens,  params)
  return(probs)

}
