gridSearchR <- function(DF, label, group, cv_groups = 5, group.column = "group", vars, distribution = c("bernoulli", "gaussian", "poisson", "huberized"), eta, max_depth, n.minobsinnode, nrounds, method = "cv", cv.folds = 4, cl = NULL) {
  model <- as.formula(paste0(label, "~",
                             paste(vars, collapse = "+")))
  GRP <- group
  set.seed(1)
  TRAIN <- lapply(1:cv_groups, function(i) {
    DF[DF[[group.column]] != i, ]
  })
  TEST <- lapply(1:cv_groups, function(i) {
    DF[DF[[group.column]] == i, ]
  })
  COMBN <- expand.grid(eta, max_depth, n.minobsinnode)
  output <- lapply(split(COMBN, row.names(COMBN)), function(m) {
    if(is.null(cl)) {
      case.gbm <- lapply(1:cv_groups, function(j) {
      gbm(data=TRAIN[[j]],
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
        gbm(data=TRAIN[[j]],
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
    best.iter <- sapply(case.gbm, gbm.perf, method = method, plot.it=F)
    if(distribution %in% c("bernoulli", "huberized")) {
      output2 <- lapply(1:cv_groups, function(y) {
        cbind(predict(case.gbm[[y]],
                      newdata = TRAIN[[y]],
                      n.trees = best.iter[y],
                      type = "response"),
              as.numeric(TRAIN[[y]][, label]))
        })
      eval_train <- sapply(output2, function(x) colAUC(x[,1],x[,2]))
      rmse_train <- sapply(1:cv_groups, function(x) {
        Metrics::rmse(actual = TRAIN[[x]][, label],
                      predicted = predict(case.gbm[[x]],
                                          newdata = TRAIN[[x]],
                                          n.trees = best.iter[x],
                                          type = "response"))
        })
      output2 <- lapply(1:cv_groups, function(y) {
        cbind(predict(case.gbm[[y]], 
                      newdata = TEST[[y]],
                      n.trees = best.iter[y],
                      type = "response"),
              as.numeric(TEST[[y]][, label]))
        })
      eval_test <- sapply(output2, function(x) colAUC(x[,1],x[,2]))
      rmse_test <- sapply(1:cv_groups, function(x) {
        Metrics::rmse(actual = TEST[[x]][, label],
                      predicted = predict(case.gbm[[x]],
                                          newdata = TEST[[x]],
                                          n.trees = best.iter[x],
                                          type = "response"))
        })
    } else {
      eval_train <- sapply(1:cv_groups, function(x) {
        1-sum((TRAIN[[x]][, label] - predict(case.gbm[[x]],
                                        newdata=TRAIN[[x]],
                                        n.trees=best.iter[x],
                                        type="response"))^2)/sum((TRAIN[[x]][, label] - mean(TRAIN[[x]][, label]))^2)
        })
      rmse_train <- sapply(1:cv_groups, function(x) {
        Metrics::rmse(actual = TRAIN[[x]][, label],
                      predicted = predict(case.gbm[[x]],
                                          newdata = TRAIN[[x]],
                                          n.trees = best.iter[x],
                                          type = "response"))
        })
      eval_test <- sapply(1:cv_groups, function(x) {
        1-sum((TEST[[x]][, label] - predict(case.gbm[[x]],
                                       newdata=TEST[[x]],
                                       n.trees =best.iter[x],
                                       type="response"))^2)/sum((TEST[[x]][, label] - mean(TEST[[x]][, label]))^2)
        })
      rmse_test <- sapply(1:cv_groups, function(x) {
        Metrics::rmse(actual = TEST[[x]][, label],
                      predicted = predict(case.gbm[[x]],
                                          newdata = TEST[[x]],
                                          n.trees = best.iter[x],
                                          type = "response"))
        })
    }
    vals <- cbind.data.frame(eta = m[1],
                       max_depth = m[2],
                       n.minobsinnode = m[3],
                       n.trees = best.iter,
                       eval_train = eval_train,
                       eval_test = eval_test,
                       rmse_train = rmse_train,
                       rmse_test = rmse_test,
                       group = 1:cv_groups)
    error_plot <- cbind.data.frame(index = 1:nrounds,
                       train_mean = colMeans(do.call(rbind, lapply(case.gbm, function(x) x$train.error))),
                       train_sd = apply(do.call(rbind, lapply(case.gbm, function(x) x$train.error)), 2, sd),
                       valid_mean = colMeans(do.call(rbind, lapply(case.gbm, function(x) x$cv.error))),
                       valid_sd = apply(do.call(rbind, lapply(case.gbm, function(x) x$cv.error)), 2, sd),
                       best.iter = round(mean(best.iter)),
                       group = paste0("eta:",
                                      m[1],
                                      ", max depth:",
                                      m[2],
                                      ", min obs in node:",
                                      m[3]))
    list(vals, error_plot)
    })
    list(do.call(rbind, lapply(output, function(i) i[[1]])),
         do.call(rbind, lapply(output, function(j) j[[2]])))
}
