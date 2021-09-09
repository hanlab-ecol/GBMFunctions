gbmCrossVal <- function(cv.folds, nTrain, n.cores,
                        class.stratify.cv, data,
                        x, y, offset, distribution, w, var.monotone,
                        n.trees, interaction.depth, n.minobsinnode,
                        shrinkage, bag.fraction,
                        var.names, response.name, group) {
  i.train <- 1:nTrain
  cv.group <- GRP#getCVgroup(distribution, class.stratify.cv, y,
  #                       i.train, cv.folds, group)
  ## build the models
  cv.models <- gbmCrossValModelBuild(cv.folds, cv.group, n.cores,
                                     i.train, x, y, offset,
                                     distribution, w, var.monotone,
                                     n.trees, interaction.depth,
                                     n.minobsinnode, shrinkage,
                                     bag.fraction, var.names,
                                     response.name, group)
  
  ## get the errors
  cv.error  <- gbmCrossValErr(cv.models, cv.folds, cv.group, nTrain, n.trees)
  best.iter.cv <- which.min(cv.error)
  
  ## get the predictions
  predictions <- gbmCrossValPredictions(cv.models, cv.folds, cv.group,
                                        best.iter.cv, distribution,
                                        data[i.train, ], y)
  list(error = cv.error, predictions = predictions)
}
