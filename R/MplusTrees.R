#' Recursive partitioning trees with Mplus models
#'
#' Generates recursive partitioning trees using M\emph{plus} models. \code{MplusTrees()} takes an
#' M\emph{plus} model written in the form of an \code{MplusAutomation} script, uses
#' \code{MplusAutomation} to fit the model in M\emph{plus}, and performs recursive partitioning
#' using \code{rpart}.
#'
#' @param script An \code{MplusAutomation} script file
#' @param data Dataset that is specified in the script
#' @param rPartFormula Formula of the form ~ variable names
#' @param catvars Vector of names of categorical covariates
#' @param group id variable. If not specified an id variable is created for each row.
#' @param control Control object for \code{rpart}.
#' @param se Whether to print standard errors and \emph{p} values. In general should be set to FALSE.
#' @param psplit Whether to use likelihood ratio \emph{p} values as a splitting criterion
#' @param palpha Type I error rate (alpha level) for rejecting with likelihood ratio test when
#' \code{psplit} set to \code{TRUE}
#' @details In the fitted object, \code{rpart_out} provides the tree structure, \code{terminal} gives a
#' vector of terminal nodes, \code{where} shows the terminal node of each id, and \code{estimates} gives
#' the parameter estimates for each terminal node.
#'
#' By default \code{MplusTrees()} only splits on the criteria specified in the \code{control}
#' argument, the most important of which is the \code{cp} parameter. However, the user can also split on the
#' \emph{p} value generated from the likelihood ratio test comparing the parent node to a multiple group
#' model consisting of 2 groups (the daughter nodes). This \emph{p} value criterion is used in addition
#' to the \code{cp} criterion in that both must be met for a split to be made. The \code{psplit} argument
#' turns this option on, \code{palpha} sets the alpha level criterion for rejection.
#' @import MplusAutomation
#' @import rpart
#' @import nlme
#' @importFrom stats terms as.formula model.matrix pchisq
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
#' }

MplusTrees <- function(script,
                       data,
                       rPartFormula,
                       catvars = NULL,
                       group= ~ id,
                       control = rpart.control(),
                       se = F,
                       psplit = F,
                       palpha = .05){

  while(grepl("model_dir",getwd())){setwd("..")}
  if(dir.exists("model_dir")){unlink("model_dir",recursive = T,force=T)}
  if(!dir.exists("model_dir")){dir.create("model_dir")}
  setwd("model_dir")

  if(missing(group)){data$id = 1:nrow(data)}

  indicator = function(data,catvars=catvars,rPartFormula=rPartFormula){
    if(any((catvars %in% names(data))==F)){
      stop("catvars variable not in dataset")
    }
    catvar.nums = which(names(data)%in%catvars)
    newdata = data
    for(i in catvar.nums){
      newdata[,i] = as.factor(data[,i])
      if(length(levels(newdata[,i]))<3){
        stop("catvars variable has fewer than 3 levels")
      }
      newdata = cbind(newdata, model.matrix(with(newdata,
                      as.formula(paste("~",names(data)[i],"-1",sep="")))))
    }
    newdata=newdata[,-catvar.nums]
    covs1 = strsplit(as.character(rPartFormula)[2]," \\+ ")[[1]]
    if(any((catvars %in% covs1)==F)){
      stop("catvars variable not in rPartFormula")
    }
    covs2 = names(newdata)[(ncol(data)-length(catvars)+1):ncol(newdata)]
    covs3 = c(covs1,covs2)
    covs4 = covs3[-which(covs3%in%catvars)]
    rpf = as.formula(paste("~",paste(covs4,collapse="+"),sep=""))

    catdata = vector("list")
    catdata$data = newdata
    catdata$rpf = rpf
    return(catdata)
  }

  if(!is.null(catvars)){
    catdata = indicator(data,catvars=catvars,rPartFormula=rPartFormula)
    data = catdata$data
    rPartFormula = catdata$rpf
  }

  groupingName = attr(terms(splitFormula(group,'~')[[1]]),"term.labels")
  groupingFactor = data[,names(data)==groupingName]


  evaluation <- function(y, wt, parms){
    script$rdata = data=parms[groupingFactor%in%y,]
    fit = mplusModeler(script,run=1L,modelout = "Model.1.inp")
    fitt=5
    list(label=fitt,deviance=-2*(fit$results$summaries$LL))
  }


  split.tree <- function(y, wt, x, parms, continuous) {
    print(paste("splitting:", length(unique(x)), "values"))
    dev = vector()
    xUnique = unique(x)
    script$rdata = data=parms[groupingFactor%in%y,]
    fit = mplusModeler(script,run=1L,modelout = "Model.1.inp")
    rootDev = fit$results$summaries$LL

    mplusmodwrap = function(index,script,run,modelout){
      if(index == 1){
        script$rdata = parms[groupingFactor %in% yLeft, ]
        out = mplusModeler(script,run=run,modelout=modelout)
      }
      if(index == 2){
        script$rdata = parms[groupingFactor %in% yRight, ]
        out = mplusModeler(script,run=run,modelout=modelout)
      }
      return(out)
    }

    if (continuous) {
      for (i in xUnique) {
        yLeft = y[x <= i]
        yRight = y[x > i]
        if (length(yLeft) < control$minbucket || length(yRight) <
            control$minbucket) {
          dev = c(dev, 0)
        }
        else{
          script$rdata = parms[groupingFactor %in% yLeft, ]
          modelLeft = mplusModeler(script,run=1L,modelout = "Model.1.inp")
          script$rdata = parms[groupingFactor %in% yRight, ]
          modelRight = mplusModeler(script,run=1L,modelout = "Model.1.inp")

          if(psplit==T){
            chisq = -2*(fit$results$summaries$LL-
              modelLeft$results$summaries$LL-modelRight$results$summaries$LL)
            p0 = fit$results$summaries$NIndependentVars +
              fit$results$summaries$NDependentVars
            df = p0^2 + 3*p0 - fit$results$summaries$Parameters -
              modelLeft$results$summaries$ChiSqM_DF -
              modelRight$results$summaries$ChiSqM_DF
            pval = 1 - pchisq(chisq,df)
            if(length(pval)!=0){
              if(pval >= palpha){
                dev = c(dev, 0)
              }
              else if(pval < palpha){
                dev = c(dev, modelLeft$results$summaries$LL +
                          modelRight$results$summaries$LL)
              }
            }
            else{
              dev = c(dev, modelLeft$results$summaries$LL +
                        modelRight$results$summaries$LL)
            }
          }
          else if(psplit==F){
            dev = c(dev, modelLeft$results$summaries$LL +
                      modelRight$results$summaries$LL)
          }
        }
      }
      good = rep(0, length(x))
      for (i in 1:length(xUnique)) {
        good[x == xUnique[i]] = dev[i]
      }
      good = good[1:(length(good) - 1)]
      list(goodness = good + abs(rootDev) * (good != 0) * 2,
           direction = rep(-1, length(good)))
    }
    else {
      order = rep(0, length(xUnique))
      dir = sort(order, index.return = TRUE)$ix
      for (i in 1:(length(dir) - 1)) {
        yLeft = y[x %in% dir[1:i]]
        yRight = y[x %in% dir[(i + 1):length(dir)]]
        if (length(yLeft) < control$minbucket ||
            length(yRight) < control$minbucket) {
          dev = c(dev, 0)
        }
        else {
          script$rdata = parms[groupingFactor %in% yLeft, ]
          modelLeft = mplusModeler(script,run=1L,modelout = "Model.1.inp")
          script$rdata = parms[groupingFactor %in% yRight, ]
          modelRight = mplusModeler(script,run=1L,modelout = "Model.1.inp")

          if(psplit==T){
            chisq = -2*(fit$results$summaries$LL-
              modelLeft$results$summaries$LL-modelRight$results$summaries$LL)
            p0 = fit$results$summaries$NIndependentVars +
              fit$results$summaries$NDependentVars
            df = p0^2 + 3*p0 - fit$results$summaries$Parameters -
              modelLeft$results$summaries$ChiSqM_DF -
              modelRight$results$summaries$ChiSqM_DF
            pval = 1 - pchisq(chisq,df)
            if(length(pval)!=0){
              if(pval >= palpha){
                dev = c(dev, 0)
              }
              else if(pval < palpha){
                dev = c(dev, modelLeft$results$summaries$LL +
                          modelRight$results$summaries$LL)
              }
            }
            else{
              dev = c(dev, modelLeft$results$summaries$LL +
                        modelRight$results$summaries$LL)
            }
          }
          else if(psplit==F){
            dev = c(dev, modelLeft$results$summaries$LL +
                      modelRight$results$summaries$LL)
          }
        }
      }
      list(goodness = dev + abs(rootDev) * (dev != 0) * 2, direction = dir)
    }
  }


  initialize <- function(y,offset,parms=0,wt){
    list(
      y=y,
      parms=parms,
      numresp=1,
      numy=1,
      summary=function(yval,dev,wt,ylevel,digits){paste("deviance (-2logLik)",
                       format(signif(dev),3),"fitt",signif(yval,2))},
      text= function(yval,dev,wt,ylevel,digits,n,use.n){
        if(!use.n){paste("m:",format(signif(yval,1)))}
        else{paste("n:",n)}
      }
    )
  }
  model <- list()
  estimates <- list()
  model.rpart = rpart(paste(groupingName,c(rPartFormula)),
                      method=list(eval=evaluation, split=split.tree,init=initialize),
                      control=control,data=data,parms=data)
  model$rpart_out <- model.rpart
  frame <- model.rpart$frame
  node2 = row.names(frame[frame[,"var"] == "<leaf>",])
  model$terminal = node2

  model$where = model.rpart$where
  oldwhere = names(table(model$where))
  owind = vector("list",length(oldwhere))
  for(i in 1:length(oldwhere)){
    owind[[i]] = which(model$where == oldwhere[i])
  }
  for(i in 1:length(oldwhere)){
    model$where = replace(model$where,owind[[i]],node2[i])
  }

  for(j in 1:length(table(model.rpart$where))){
    id <- names(table(model.rpart$where))[j]==model.rpart$where
    script$rdata = data[id,]
    fit = mplusModeler(script,run=1L,modelout = "Model.1.inp")
    if(se==F){estimates[[node2[j]]] = fit$results$parameters[[1]][,1:3]}
    else if(se==T){estimates[[node2[j]]] = fit$results$parameters}
  }
  model$estimates = estimates

  class(model) <- "mplustree"

  setwd("..")

  return(model)
}
