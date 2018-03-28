#' Mplus model recursive partitioning trees
#'
#' @param script An MplusAutomation script file
#' @param data Dataset that is specified in the script
#' @param rPartFormula Formula of the form ~ variable names
#' @param group What variable is an id variable. If not specified
#'              an id variable is created for each row.
#' @param parallel Whether to parallelize each node (max two).
#' @param control Control object for rpart.
#' @import MplusAutomation
#' @import rpart
#' @import snowfall
#' @import nlme
#' @importFrom stats terms
#' @export


MplusTrees <- function(script,
                data,
                rPartFormula,
                group= ~ id,
                parallel=FALSE,
                control = rpart.control()){

  # have to have id variable
  dir.create("model_dir")
  setwd("model_dir")

  MplusAutomation=NULL

  if(parallel==TRUE){
    sfInit(T,2)
    sfLibrary(package=MplusAutomation)
  }

  groupingName = attr(terms(splitFormula(group,'~')[[1]]),"term.labels")
  groupingFactor = data[,names(data)==groupingName]
#    terms = attr(terms(lmeFormula),"term.labels")
   # continuous = !is.factor(data[,names(data)==terms[1]])
    ### The 3 subfunctions necessary for rpart to work.
    # The evaluation function.
    # Called once per node:
    # returns a list of two variables: a label for the node
    # and a deviance value for the node.  The deviance is
    # of length one, equal to 0 if the node is perfect/ unsplittable
    # larger equals worse
    evaluation <- function(y, wt, parms){
      script$rdata = data=parms[groupingFactor%in%y,]
      #fit = mplusModeler(script,run=1L,modelout = "Model.1.inp")
      fit = try(mplusModeler(script,run=1L,modelout = "Model.1.inp"),silent=F)
      fitt=5
      list(label=fitt,deviance=-2*(fit$results$summaries$LL))
    }


  # The split function, where the work occurs.  This is used to decide
  # on the optimal split for a given covariate.
  # Called once per split candidate per node
  ### If the covariate, x, is continuous:
  # x variable is ordered
  # y is provided in the sort order of x
  # returns two vectors of length (n-1)
  #      goodness: goodness of the split, with larger numbers better. 0=no split
  #      direction: -1 = send y < cutpoint to the left
  #                  1 = send y < cutpoint to the right
  #
  ### If x is non-continuous
  # x is a set of integers (NOT the original values) defining the
  # groups for an unordered predictor.
  # Again, return goodness and direction
  #      direction: a vector of length m (number of groups), which is the applied
  #                 ordering for the unordered categories.  This is done so that
  #                 m-1 splits are performed, instead of all possible splits
  #      goodness: m-1 values, same idea as before

  ### pass in the dataset through the parms variable, with subj as y
  split <- function(y, wt, x, parms, continuous) {
    print(paste("splitting:", length(unique(x)), "values"))
    dev = vector()
    xUnique = unique(x)

    #  rootDev = lme(lmeFormula, data = parms[groupingFactor %in%
    #                                           y, ], random = randomFormula, correlation = R,
    #                na.action = na.omit)$logLik
    script$rdata = data=parms[groupingFactor%in%y,]
     #fit = mplusModeler(script,run=1L,modelout = "Model.1.inp")
     fit = try(mplusModeler(script,run=1L,modelout = "Model.1.inp"),silent=F)
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
        else {

          if(parallel==TRUE){
            modelLR = sfClusterApplyLB(1:2, mplusmodwrap, script,run=1L,modelout = "Model.1.inp")

            modelLeft = modelLR[[1]]
            modelRight = modelLR[[2]]
          }else{
            script$rdata = parms[groupingFactor %in%
                                   yLeft, ]
            modelLeft = try(mplusModeler(script,run=1L,modelout = "Model.1.inp"),silent=F)
            if(inherits(modelLeft, "try-error")) {
              modelLeft = NA
            }else{
              modelLeft = modelLeft
            }
            script$rdata = parms[groupingFactor %in%
                                   yRight, ]
            modelRight = try(mplusModeler(script,run=1L,modelout = "Model.1.inp"),silent=F)
            if(inherits(modelRight, "try-error")) {
              modelRight = NA
            }else{
              modelRight = modelRight
            }
          }


          if(is.null(modelLeft$results$summaries$LL) == T|
             is.null(modelRight$results$summaries$LL) == T){
            dev = c(dev,-9e10)
          }else{
            dev = c(dev, modelLeft$results$summaries$LL + modelRight$results$summaries$LL)
          }



        }
      }
      good = rep(0, length(x))
      for (i in 1:length(xUnique)) {
        good[x == xUnique[i]] = dev[i]
      }
      good = good[1:(length(good) - 1)]
      list(goodness = good + abs(rootDev) * (good != 0) *
             2, direction = rep(-1, length(good)))
    }
    else {
      order = rep(0, length(xUnique))
     # response = parms[, names(parms) == responseName]  #  -          ----- problem
     # for (i in 1:length(xUnique)) {
     #   order[i] = mean(response[x == xUnique[i]], na.rm = TRUE)
     # }
      dir = sort(order, index.return = TRUE)$ix
      for (i in 1:(length(dir) - 1)) {
        yLeft = y[x %in% dir[1:i]]
        yRight = y[x %in% dir[(i + 1):length(dir)]]
        if (length(yLeft) < control$minbucket || length(yRight) <
            control$minbucket) {
          dev = c(dev, 0)
        }
        else {


          if(parallel==TRUE){
           # sfInit(T,2)
           # sfLibrary(package=MplusAutomation)
            modelLR = sfClusterApplyLB(1:2, mplusmodwrap, script,run=1L,modelout = "Model.1.inp")

            modelLeft = modelLR[[1]]
            modelRight = modelLR[[2]]
          }else{
            script$rdata = parms[groupingFactor %in%
                                   yLeft, ]
            modelLeft = try(mplusModeler(script,run=1L,modelout = "Model.1.inp"),silent=F)
            if(inherits(modelLeft, "try-error")) {
              modelLeft = NA
            }else{
              modelLeft = modelLeft
            }
            script$rdata = parms[groupingFactor %in%
                                   yRight, ]
            modelRight = try(mplusModeler(script,run=1L,modelout = "Model.1.inp"),silent=F)
            if(inherits(modelRight, "try-error")) {
              modelRight = NA
            }else{
              modelRight = modelRight
            }

          }

          if(is.null(modelLeft$results$summaries$LL) == T|
             is.null(modelRight$results$summaries$LL) == T){
            dev = c(dev,-9e10)
          }else{
            dev = c(dev, modelLeft$results$summaries$LL + modelRight$results$summaries$LL)
          }

        }
      }
      list(goodness = dev + abs(rootDev) * (dev != 0) *
             2, direction = dir)
    }
  }
  # The init function.  This is used, to the best of my knowledge, to initialize the process.
  # summary is used to fill print the report summary(model), and text is used to add text to
  # the plot of the tree.
  initialize <- function(y,offset,parms=0,wt){
    list(
      y=y,
      parms=parms,
      numresp=1,
      numy=1,
      summary=function(yval,dev,wt,ylevel,digits){paste("deviance (-2logLik)",format(signif(dev),3),"fitt",signif(yval,2))},
      text= function(yval,dev,wt,ylevel,digits,n,use.n){
        if(!use.n){paste("m:",format(signif(yval,1)))}
        else{paste("n:",n)}
      }
    )
  }
  model <- list()
  summary <- list()
  model.rpart = rpart(paste(groupingName,c(rPartFormula)),
                      method=list(eval=evaluation,
                      split=split,init=initialize),
                      control=control,data=data,parms=data)

  model$rpart_out <- model.rpart

  frame <- model.rpart$frame
  node2 = row.names(frame[frame[,"var"] == "<leaf>",])
  model$node2 = node2

  for(j in 1:length(table(model.rpart$where))){
    id <- names(table(model.rpart$where))[j]==model.rpart$where

    script$rdata = data[id,]
    fit = mplusModeler(script,run=1L,modelout = "Model.1.inp")

    summary[[node2[j]]] = fit$results$parameters


  }
if(parallel==TRUE){
  sfStop()
}
  model$summary = summary
 unlink(getwd(),recursive=T)
  model
}
