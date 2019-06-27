#' Plots tree structure of an Mplus Tree
#'
#' Wrapper using \code{rpart.plot} package to plot the tree structure of a
#' fitted Mplus Tree
#'
#' @param x An object of class "mplustree" (a fitted Mplus Tree)
#' @param ... Other arguments passed to \code{rpart.plot}
#' @method plot mplustree
#' @import rpart.plot
#' @details Each node of the plot by default contain the -2
#'          log-likelihood (deviance), the number of individuals
#'          in the node, and the percentage of the total sample
#'          in the node.
#' @author Sarfaraz Serang, relying heavily on the \code{rpart.plot}
#' package by Stephen Milborrow.
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
#' fit
#'
#' plot(fit)
#' }

plot.mplustree <- function(x, ...){
  if(is.null(list(...)$extra)){extra = 101}
  if(is.null(list(...)$box.pallete)){box.palette=0}

  x$rpart_out[[1]]$yval=x$rpart_out[[1]]$dev

  suppressWarnings(rpart.plot::rpart.plot(
    x$rpart_out,extra=extra,
    box.palette=box.palette,...))
}
