#' Summarizing MplusTrees model Fits
#'
#' \code{summary} method for class "mplustree".
#'
#' @param object An object of class "mplustree" (a fitted Mplus Tree)
#' @param ... Other arguments passed to or from other methods
#' @method summary mplustree
#' @details Prints the tree structure given in \code{object}
#' @export
#' @examples
#' \dontrun{
#' library(lavaan)
#'
#' script = mplusObject(
#'    TITLE = "Example #1 - Factor Model;",
#'    MODEL = "f1 BY x1-x3; f2 BY x4-x6; f3 BY x7-x9;",
#'    usevariables = c('x1','x2','x3','x4','x5','x6','x7','x8','x9'),
#'    rdata = HolzingerSwineford1939)
#'
#' fit = MplusTrees(script, HolzingerSwineford1939, group=~id,
#'    rPartFormula=~sex+school+grade, control=rpart.control(cp=.01))
#'
#' summary(fit)
#' }

summary.mplustree <- function(object, ...){
  return(object$rpart_out)
}
