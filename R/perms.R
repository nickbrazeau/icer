#' Function to create all moi permutations
#'
#' @param m Total number we are summing to
#' @param n Number of different non-negative integers used. Default = 3
#'
#' @details Creates a matrix of the different ways of adding
#'   to m with n non-negative integers.
#'
#' @return Returns a matrix with n columns, with each row summing to m.
#'
#' @examples
#' \dontrun{
#'
#' moi_perms(2,3)
#'
#' icer:::moi_perms(m = 2, n = 3)
#'
#' #>      [,1] [,2] [,3]
#' #> [1,]    0    0    2
#' #> [2,]    1    0    1
#' #> [3,]    2    0    0
#' #> [4,]    0    1    1
#' #> [5,]    1    1    0
#' #> [6,]    0    2    0
#' }

moi_perms <- function(m, n = 3) {
  expr <- parse(text = paste0("expand.grid(", paste(rep("0:m",
                                                        n - 1), collapse = ","), ")"))
  if(n-1 == 1){
    X <- t(as.matrix(eval(expr)))
    X <- t( rbind(X, X[rev(1:ncol(X))]) )
    colnames(X) <- NULL
    return(X)
  } else{
    X <- t(as.matrix(eval(expr)))
    X <- X[, colSums(X) <= m]
    X <- t(rbind(X, m - colSums(X)))
    colnames(X) <- NULL
    return(X)
  }

}

#' Determines the potential compositions from the cases names
#' (e.g. composition levels)
#'
#' @param names character vector of case names
#'
#' @details Creates a named vector of all combinations for the names
#'   in probs, separated by a "/" and sorted by character length
#'   and then alphabetically
#'
#' @return Returns a vector of composition levels
#'
#' @examples
#' \dontrun{
#'
#' obj <- new("casedetect")
#' obj@cases <- c(50, 50, 100)
#' obj@casenames <- c("a", "b", "c")
#'
#' icer:::generate_composition_levels(names = obj@casenames)
#'
#' #> a     b     c   a/b   a/c   b/c a/b/c
#'
#' }

generate_composition_levels <- function(names) {

  # first lets create all the possible combinations
  levels <- parse(text = paste0("expand.grid(",
                                paste0(rep(
                                  "c(TRUE,FALSE)", length(names)
                                ), collapse = ","),
                                ")"))
  levels <- head(eval(levels), 2 ^ length(names) - 1)

  # tidy the names of the levels
  levels <-
    (apply(levels, 1, function(x) {
      paste0(names[x], collapse = "/")
    }))
  levels <- levels[order(nchar(levels), levels)]

  # return composition names
  return(levels)

}



#' Creates composition levels from the infection compositions
#'
#' @param nms Names of the species we are looking at.
#' @param perms List of outputs of \code{moi_perms}.
#'
#' @details For each matrix of moi permutations within `perms`, each row
#'   is used to work out the infection type from the infection composition.
#'   E.g. if the row of the matrix is [2, 1, 0] and our `nms` are [a, b, c]
#'   then the infeciton composition is "a/a/b" but the observed infection
#'   type will be "a/b"
#'
#' @return Returns a list of vectors that correspond to the infection types
#'   of the infection compostions given within `perms`
#'
#' @examples
#' \dontrun{
#'
#'
#' nms <- obj@casenames
#' perms <- lapply(1:2, moi_perms, n = 3)
#'
#' find_composition_level(nms,perms)
#'
#' #> [[1]]
#' #> [1] "c" "a" "b"
#' #>
#' #> [[2]]
#' #> [1] "c"   "a/c" "a"   "b/c" "a/b" "b"
#'
#' }

find_composition_level <- function(nms, perms) {

  infs <- lapply(perms, function(x) {
    apply(x, 1, function(x) { paste0(nms[x > 0], collapse = "/") })
  })

  return(infs)

}



