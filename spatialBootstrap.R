bootstrapGBMR <- function(DF, label, group, cv_groups = 5, group.column = "group", vars, k_split, distribution = c("bernoulli", "gaussian", "poisson", "huberized"), eta, max_depth, n.minobsinnode, nrounds, nruns, bootstrap = c("observed", "null"), method = "cv", cv.folds = 4, cl = NULL) {
  model<-as.formula(paste(label, "~",
                          paste(vars, collapse = "+"),
                          sep = ""))
  TRAIN <- lapply(1:cv_groups, function(i) {
    DF[DF[[group.column]] != i, ]
  })
  if(bootstrap == "null") {
    for(i in 1:cv_groups) {
      TRAIN[[i]][, label] <- sample(x = TRAIN[[i]][, label], size = nrow(TRAIN[[i]]), replace = F)
    }
  }
  TEST <- lapply(1:cv_groups, function(i) {
    DF[DF[[group.column]] == i, ]
  })
  if(bootstrap == "null") {
    for(i in 1:cv_groups) {
      TEST[[i]][, label] <- sample(x = TEST[[i]][, label], size = nrow(TEST[[i]]), replace = F)
    }
  }
  if(is.null(cl)) {
    case.gbm <- lapply(1:cv_groups, function(m) {
      gbm(data= TRAIN[[m]],
          model,
          distribution = distribution,
          n.trees = nrounds,
          shrinkage = eta,
          cv.folds = cv.folds,
          interaction.depth = max_depth,
          n.minobsinnode = n.minobsinnode,
          bag.fraction = 0.5,
          verbose = FALSE)
    })
  } else {
    clusterExport(cl, list("gbm"))
    case.gbm <- clusterApply(cl, 1:cv_groups, function(m) {
      gbm(data=TRAIN[[m]],
          model,
          distribution = distribution,#default
          n.trees = nrounds,
          shrinkage = eta,
          cv.folds = cv.folds,
          interaction.depth = max_depth,
          n.minobsinnode = n.minobsinnode,
          bag.fraction = 0.5,
          verbose = FALSE)
    })
  }
best.iter <- sapply(case.gbm, gbm.perf, method = method, plot.it=F) #this gives you the optimal number of trees 
  for(i in 1:cv_groups) {
    case.gbm[[i]]$var.levels <- lapply(case.gbm[[i]]$var.levels, function(x) replace(x, is.infinite(x), 0))
  }
  predictions <- do.call(rbind, lapply(1:cv_groups, function(j) {
    data.frame(predictions = predict(case.gbm[[j]],
                                     newdata = DF,
                                     n.trees = best.iter[j],
                                     type = "response"),
               group_run = j,
               original_value = DF[, label])
    }))
  if(distribution %in% c("bernoulli", "huberized")) {
    output2 <- lapply(1:cv_groups, function(x) cbind(predict(case.gbm[[x]],
                                                                    newdata = TRAIN[[x]],
                                                                    n.trees = best.iter[x],
                                                                    type = "response"),
                                                            as.numeric(TRAIN[[x]][, label])))
    eval_train <- sapply(output2, function(x) colAUC(x[,1],x[,2]))
    rmse_train <- sapply(1:cv_groups, function(x) {
      Metrics::rmse(actual = TRAIN[[x]][, label],
                    predicted = predict(case.gbm[[x]],
                                        newdata = TRAIN[[x]],
                                        n.trees = best.iter[x],
                                        type = "response"))
      })
    output2 <- lapply(1:length(case.gbm), function(x) {
      cbind(predict(case.gbm[[x]],
                    newdata = TEST[[x]],
                    n.trees = best.iter[x],
                    type = "response"),
            as.numeric(TEST[[x]][, label]))
      })
    eval_test <- sapply(output2, function(x) colAUC(x[,1],x[,2]))
    rmse_test <- sapply(1:length(case.gbm), function(x) {
      Metrics::rmse(actual = TEST[[x]][, label],
                    predicted = predict(case.gbm[[x]],
                                        newdata = TEST[[x]],
                                        n.trees = best.iter[x],
                                        type = "response"))
      })
  } else {
    eval_train <- sapply(1:length(case.gbm), function(x) {
      1-sum((TRAIN[[x]][, label] - predict(case.gbm[[x]],
                                               newdata=TRAIN[[x]],
                                               n.trees=best.iter[x],
                                               type="response"))^2)/sum((TRAIN[[x]][, label] - mean(TRAIN[[x]][, label]))^2)
      })
    rmse_train <- sapply(1:length(case.gbm), function(x) {
      Metrics::rmse(actual = TRAIN[[x]][, label],
                    predicted = predict(case.gbm[[x]],
                                        newdata = TRAIN[[x]],
                                        n.trees = best.iter[x],
                                        type = "response"))
      })
    eval_test <- sapply(1:length(case.gbm), function(x) {
      1-sum((TEST[[x]][, label] - predict(case.gbm[[x]],
                                               newdata=TEST[[x]],
                                               n.trees =best.iter[x],
                                               type="response"))^2)/sum((TEST[[x]][, label] - mean(TEST[[x]][, label]))^2)
      })
    rmse_test <- sapply(1:length(case.gbm), function(x) {
      Metrics::rmse(actual = TEST[[x]][, label],
                    predicted = predict(case.gbm[[x]],
                                        newdata = TEST[[x]],
                                        n.trees = best.iter[x],
                                        type = "response"))
      })
  }
  df_importance <- do.call(rbind, lapply(case.gbm, function(j) t(as.data.frame(summary(j)$rel.inf[match(vars, summary(j)$var)], optional = T))))
  colnames(df_importance) <- vars
  pd_out <- foreach(i = 1:cv_groups, .combine = "rbind") %do% {
    pd_out <- lapply(vars, function(m) dplyr::mutate(plot.gbm(case.gbm[[i]],
                                                              i.var = m,
                                                              return.grid = T,
                                                              type = "response"),
                                                     variable.name = m,
                                                     effect = "marginal.effect",
                                                     group_run = i))
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
                           group_run = 1:cv_groups,
                           df_importance)
  list(out1, pd_out, predictions)
}
