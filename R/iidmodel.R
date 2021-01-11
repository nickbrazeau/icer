#' Independent model of strain acquisition
#'
#' @param obj S4 object;
#' @param params Vector of parameters including the frequency of each species
#'   and the mean moi mu.
#' @param data Real data to compare model predictions to.
#' @param pci_list Permutation, composition and infections list.
#' @param max_moi Maxuimum infction composition explored. Default = 25.
#' @details For a given frequency of each strain and mean moi, returns the
#'   related multionomial distribution describing each infection type.
#'
#' @return Returns a named multinomial probability distribution.
independent <- function(obj, params, data, pci_list,  max_moi = 25) {

  # handle frequency params
  # ------------------------------------
  n_sp <- ncol(pci_list$perms[[1]])
  freqs <- params[seq_len(n_sp)]
  freqs <- freqs / sum(freqs)

  # Density of each infection composition using a multinomial distribution
  densities <- lapply(pci_list$perms, function(x){
    dens <- apply(x, 1, dmultinom, size = sum(x[1,]), prob = freqs)
  })

  # turn these into a a matrix of multinomial densities
  multi_dens <- comp_dens_to_multi_dens(pci_list = pci_list, densities = densities)
  w <- moi_probabilities(obj = obj, maxmoi = max_moi, params = params)

  # calculate the overall multinomial distributon
  probs <- get_join_comp_moi_prob(obj = obj, multi_dens = multi_dens, w = w)

  return(probs)

}



