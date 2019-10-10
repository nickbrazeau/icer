#' Survelliance Class
#'
#' @details Class for surveillance data, or when the the true prevalence
#' is known. For a true prevalence to be known, we must have the
#' population at risk as the denominator (e.g. all individuals
#' that could have been sampled).
#'

setClass("surveillance", slots = list(
  denominator = "numeric",
  cases = "numeric",
  casenames = "character"

))


#' Case-Detection Class
#'
#' @details Class for case-detection data where the true
#' prevalence of cases is not known. Instead, we just have the
#' number of cases that were recorded.

setClass("casedetect", slots = list(
  cases = "numeric",
  casenames = "character"
))




#' Check to see if classes were made correctly
#'
#' @param obj S4 object;
#' @return logical for whether class was created correctly or not
#' @noRd
#'

checkclass <- function(obj){
  # note, no assertions needed to check slot types
  # this is handled by S4 create
  assert_custom_class(x = obj, c = c("surveillance", "casedetect"))
  if(assert_custom_class(x = obj, c = c("surveillance"))){

    assert_length(obj@denominator, 1)
    assert_gr(x = obj@denominator, y = 0)
    assert_gr(x =sum(obj@cases), y = 0)
    assert_leq(x = sum(obj@cases), y = obj@denominator)
    assert_same_length(obj@cases, obj@casenames)


  } else if(assert_custom_class(x = obj, c = c("casedetect"))){
    assert_same_length(obj@cases, obj@casenames)
    assert_gr(x =sum(obj@cases), y = 0)
  }

}

