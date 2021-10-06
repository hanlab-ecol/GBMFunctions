######Bootstrap function for gbm
##The bootstrapGBM function takes a variety of variables:
#DF - a data.frame of your label and predictor variables
#label - column name of your label column
#vars - predictor variable column names (may be easiest to use the colnames function and index out the label column)
#k_split -  percentage of the data to use for training
#distribution to use for gbm (currently tested with bernoulli, gaussian, and poisson), bernoulli is the default
#eta - learning rate determined from parameterization to use in bootstrapping
#max_depth - depth rate determined from parameterization
#n.minobsinnode - vector of numbers of minimum observations in node to iterate through
#nrounds - max number of trees to set
#nruns - number of bootstrap runs to use
#bootstrap - determines whether to use the observed data or a null distribution for performing the bootstrap
#method - set at cv
#cv.folds - number of folds to use with gbm
#cl - default is NULL. This is a numeric vector of length 2 with the first being the number of cores used to run gbm 
#within clusterApply and the second supplying the n.cores argument for gbm::gbm
######
##Example
#packages <- c("gbm", "caret", "Matrix", "pdp", "caTools", "ROCR", "dplyr", "foreach", "dismo", "doSNOW", "parallel")
#sapply(packages, library, character.only = T)
#detectCores()
#cl <- makeCluster(20, "SOCK")
#registerDoSNOW(cl)
#alt <- bootstrap_gbm(DF, label = "sr", vars = colnames(DF)[-1:-3], eta = 0.001, max_depth = 2, nrounds = 1000, distribution = "gaussian", k = 5, nruns = 5, bootstrap = "observed")
######
bootstrapGBM <- function(DF, label, vars, k_split, distribution = c("bernoulli", "gaussian", "poisson", "huberized"), eta, max_depth, n.minobsinnode, nrounds, nruns, bootstrap = c("observed", "null"), method = "cv", cv.folds = 5, cl = NULL) {
  model<-as.formula(paste(label, "~",
                          paste(vars, collapse = "+"),
                          sep = ""))
  DaFr <- foreach(i = 1:nruns) %do% {
    set.seed(i)
    if(distribution == "bernoulli") {
      DP <- createDataPartition(as.factor(DF[, label]), p = k_split)[[1]]
    } else DP <- createDataPartition(DF[, label], p = k_split)[[1]]
    if(bootstrap == "observed") {
      TRAIN <- DF[DP, ]
      TEST <- DF[-DP, ]
    } else {#"null"
      TRAIN <- DF[DP, ]
      TRAIN[, label] <- sample(x = TRAIN[, label], size = nrow(TRAIN), replace = F)
      TEST <- DF[-DP, ]
      TEST[, label] <- sample(x = TEST[, label], size = nrow(TEST), replace = F)
    }
    list(TRAIN, TEST)
  }
  if(is.null(cl)) {
    case.gbm <- lapply(1:nruns, function(m) {
      gbm(data= DaFr[[m]][[1]],
          model,
          distribution = distribution,
          n.trees = nrounds,
          shrinkage = eta,
          cv.folds = cv.folds,
          interaction.depth = max_depth,
          n.minobsinnode = n.minobsinnode,
          bag.fraction = 0.5,
          verbose = FALSE,
          n.cores = 1)
    })
  } else {
    cores <- makeCluster(cl[1])
    on.exit(parallel::stopCluster(cores))
    clusterExport(cores, list("gbm"))
    case.gbm <- clusterApply(cores, 1:nruns, function(m) {
      gbm(data= DaFr[[m]][[1]],
          model,
          distribution = distribution,
          n.trees = nrounds,
          shrinkage = eta,
          cv.folds = cv.folds,
          interaction.depth = max_depth,
          n.minobsinnode = n.minobsinnode,
          bag.fraction = 0.5,
          verbose = FALSE,
          n.cores = cl[2])
    })
  }
  best.iter <- sapply(case.gbm, gbm.perf, method = method, plot.it=F) #this gives you the optimal number of trees 
  for(i in 1:nruns) {
  case.gbm[[i]]$var.levels <- lapply(case.gbm[[i]]$var.levels, function(x) replace(x, is.infinite(x), 0))
  }
  predictions <- do.call(rbind, lapply(1:nruns, function(j) {
    data.frame(predictions = predict(case.gbm[[j]], newdata = DF, n.trees = best.iter[j], type = "response"),
               bootstrap_run = j,
               original_value = DF[, label])}))
  if(distribution %in% c("bernoulli", "huberized")) {
    output2 <- lapply(1:length(case.gbm), function(x) cbind(predict(case.gbm[[x]],
                                                                    newdata = DaFr[[x]][[1]],
                                                                    n.trees = best.iter[x],
                                                                    type = "response"),
                                                            as.numeric(DaFr[[x]][[1]][, label])))
    eval_train <- sapply(output2, function(x) colAUC(x[,1],x[,2]))
    rmse_train <- sapply(1:length(case.gbm), function(x) Metrics::rmse(actual = DaFr[[x]][[1]][, label],
                                                                       predicted = predict(case.gbm[[x]],
                                                                                           newdata = DaFr[[x]][[1]],
                                                                                           n.trees = best.iter[x],
                                                                                           type = "response")))
    output2 <- lapply(1:length(case.gbm), function(x) cbind(predict(case.gbm[[x]],
                                                                    newdata = DaFr[[x]][[2]],
                                                                    n.trees = best.iter[x],
                                                                    type = "response"),
                                                            as.numeric(DaFr[[x]][[2]][, label])))
    eval_test <- sapply(output2, function(x) colAUC(x[,1],x[,2]))
    rmse_test <- sapply(1:length(case.gbm), function(x) Metrics::rmse(actual = DaFr[[x]][[2]][, label],
                                                                      predicted = predict(case.gbm[[x]],
                                                                                          newdata = DaFr[[x]][[2]],
                                                                                          n.trees = best.iter[x],
                                                                                          type = "response")))
  } else {
    eval_train <- sapply(1:length(case.gbm), function(x) 1-sum((DaFr[[x]][[1]][, label] - predict(case.gbm[[x]],
                                                                                                  newdata=DaFr[[x]][[1]],
                                                                                                  n.trees=best.iter[x],
                                                                                                  type="response"))^2)/sum((DaFr[[x]][[1]][, label] - mean(DaFr[[x]][[1]][, label]))^2))
    rmse_train <- sapply(1:length(case.gbm), function(x) Metrics::rmse(actual = DaFr[[x]][[1]][, label],
                                                                       predicted = predict(case.gbm[[x]],
                                                                                           newdata = DaFr[[x]][[1]],
                                                                                           n.trees = best.iter[x],
                                                                                           type = "response")))
    eval_test <- sapply(1:length(case.gbm), function(x) 1-sum((DaFr[[x]][[2]][, label] - predict(case.gbm[[x]],
                                                                                                 newdata=DaFr[[x]][[2]],
                                                                                                 n.trees =best.iter[x],
                                                                                                 type="response"))^2)/sum((DaFr[[x]][[2]][, label] - mean(DaFr[[x]][[2]][, label]))^2))
    rmse_test <- sapply(1:length(case.gbm), function(x) Metrics::rmse(actual = DaFr[[x]][[2]][, label],
                                                                      predicted = predict(case.gbm[[x]],
                                                                                          newdata = DaFr[[x]][[2]],
                                                                                          n.trees = best.iter[x],
                                                                                          type = "response")))
  }
    df_importance <- do.call(rbind, lapply(case.gbm, function(j) t(as.data.frame(summary(j)$rel.inf[match(vars, summary(j)$var)], optional = T))))
    colnames(df_importance) <- vars
    pd_out <- foreach(i = 1:nruns, .combine = "rbind") %do% {
    pd_out <- lapply(vars, function(m) dplyr::mutate(plot.gbm(case.gbm[[i]],
                                                              i.var = m,
                                                              return.grid = T,
                                                              type = "response"),
                                                     variable.name = m,
                                                     effect = "marginal.effect",
                                                     bootstrap_run = i))
    pd_out <- do.call(rbind, lapply(pd_out, function(m) {
      colnames(m)[1:2] <- c("x", "yhat")
      m}))
    pd_out}
    out1 <- cbind.data.frame(eta = eta,
                             max_depth = max_depth,
                             n.minobsinnode = n.minobsinnode,
                             n.trees = best.iter,
                             eval_train = eval_train,
                             eval_test = eval_test,
                             rmse_train = rmse_train,
                             rmse_test = rmse_test,
                             bootstrap_run = 1:nruns,
                             df_importance)
    list(out1, pd_out, predictions)
  }
