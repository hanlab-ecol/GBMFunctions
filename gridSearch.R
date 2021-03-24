######Parameterization function for gbm
##The gridSearch function takes a variety of variables:
#DF - a data.frame of your label and predictor variables
#label - column name of your label column
#vars - predictor variable column names (may be easiest to use the colnames function and index out the label column)
#k_split - percentage of the data to use for training
#distribution to use for gbm (currently tested with bernoulli, gaussian, and poisson), bernoulli is the default
#eta - learning rate vector to iterate through
#max_depth - depth vector to iterate through
#n.minobsinnode - vector of numbers of minimum observations in node to iterate through
#nrounds - max number of trees to set
#method - set at cv by default
#cv.folds - number of folds to use with gbm if your method is "cv"
#cl - the cluster object to use with parallelization. If you don't wish to use parallel processing, then this argument should be "cl = NULL"
######
##Example
#packages <- c("gbm", "caret", "Matrix", "pdp", "caTools", "ROCR", "dplyr", "foreach", "dismo", "doSNOW", "parallel")
#sapply(packages, library, character.only = T)
#detectCores()
#cl <- makeCluster(3, "SOCK")
#registerDoSNOW(cl)
#GRID <- gridSearch(TRAIT, label = "above_infection", eta = c(0.0001, 0.005, 0.001), max_depth = c(1, 2, 3, 4), n.minobsinnode = c(2, 5), vars = colnames(TRAIT)[c(8:43)], k_split = 0.8, distribution = "bernoulli", nrounds = 100000, cl = cl)
######
gridSearch <- function(DF, label, vars, k_split, distribution = c("bernoulli", "gaussian", "poisson", "huberized"), eta, max_depth, n.minobsinnode, nrounds, method = "cv", cv.folds = 5, cl) {
  model <- as.formula(paste0(label, "~",
                             paste(vars, collapse = "+")))
  set.seed(1)
  if(distribution == "bernoulli") {
    DP <- createDataPartition(as.factor(DF[, label]), p = k_split)[[1]]
  } else DP <- createDataPartition(DF[, label], p = k_split)[[1]]
  TRAIN <- DF[DP, ]
  TEST <- DF[-DP, ]
  COMBN <- expand.grid(eta, max_depth, n.minobsinnode)
  if(is.null(cl)) {
    case.gbm <- lapply(split(COMBN, row.names(COMBN)), function(m) {
      gbm(data=TRAIN,
          model,
          distribution = distribution,#default
          n.trees = nrounds,
          shrinkage = m[1],
          cv.folds = cv.folds,
          interaction.depth = m[2],
          n.minobsinnode = m[3],
          bag.fraction = 0.5,
          verbose = FALSE)
    })
  } else {
    clusterExport(cl, list("gbm"))
    case.gbm <- clusterApply(cl, split(COMBN, row.names(COMBN)), function(m) {
      gbm(data=TRAIN,
          model,
          distribution = distribution,#default
          n.trees = nrounds,
          shrinkage = m[1],
          cv.folds = cv.folds,
          interaction.depth = m[2],
          n.minobsinnode = m[3],
          bag.fraction = 0.5,
          verbose = FALSE)
    })
  }
  best.iter <- sapply(case.gbm, gbm.perf, method = method, plot.it=F) #this gives you the optimal number of trees 
  ## predictions on the TRAINING SET
  if(distribution %in% c("bernoulli", "huberized")) {
    output2 <- lapply(1:length(case.gbm), function(x) cbind(predict(case.gbm[[x]],
                                                                    newdata = TRAIN,
                                                                    n.trees = best.iter[x],
                                                                    type = "response"),
                                                            as.numeric(TRAIN[, label])))
    eval_train <- sapply(output2, function(x) colAUC(x[,1],x[,2]))
    rmse_train <- sapply(1:length(case.gbm), function(x) Metrics::rmse(actual = TRAIN[, label], predicted = predict(case.gbm[[x]], newdata = TRAIN, n.trees = best.iter[x], type = "response")))
    output2 <- lapply(1:length(case.gbm), function(x) cbind(predict(case.gbm[[x]], newdata = TEST, n.trees = best.iter[x], type = "response"), as.numeric(TEST[, label])))
    eval_test <- sapply(output2, function(x) colAUC(x[,1],x[,2]))
    rmse_test <- sapply(1:length(case.gbm), function(x) Metrics::rmse(actual = TEST[, label], predicted = predict(case.gbm[[x]], newdata = TEST, n.trees = best.iter[x], type = "response")))
  } else {
    eval_train <- sapply(1:length(case.gbm), function(x) 1-sum((TRAIN[, label] - predict(case.gbm[[x]], newdata=TRAIN, n.trees=best.iter[x], type="response"))^2)/sum((TRAIN[, label] - mean(TRAIN[, label]))^2))
    rmse_train <- sapply(1:length(case.gbm), function(x) Metrics::rmse(actual = TRAIN[, label], predicted = predict(case.gbm[[x]], newdata = TRAIN, n.trees = best.iter[x], type = "response")))
    eval_test <- sapply(1:length(case.gbm), function(x) 1-sum((TEST[, label] - predict(case.gbm[[x]], newdata=TEST, n.trees =best.iter[x], type="response"))^2)/sum((TEST[, label] - mean(TEST[, label]))^2))
    rmse_test <- sapply(1:length(case.gbm), function(x) Metrics::rmse(actual = TEST[, label], predicted = predict(case.gbm[[x]], newdata = TEST, n.trees = best.iter[x], type = "response")))
  }
  list(cbind.data.frame(eta = COMBN[, 1],
                        max_depth = COMBN[, 2],
                        n.minobsinnode = COMBN[, 3],
                        n.trees = best.iter,
                        eval_train = eval_train,
                        eval_test = eval_test,
                        rmse_train = rmse_train,
                        rmse_test = rmse_test),
       do.call(rbind, lapply(1:length(case.gbm), function(x) cbind.data.frame(index = 1:length(case.gbm[[x]]$train.error), train = case.gbm[[x]]$train.error, valid = case.gbm[[x]]$cv.error, best.iter = best.iter[x], group = paste("eta:", COMBN[x, 1], ", ", "max depth:", COMBN[x, 2], ", ", "min obs in node:", COMBN[x, 3], sep = "")))))
}
