# GBMFunctions
A collection of functions designed to help in the parameterization, evaluation, and exploration of models using {gbm}.

## Current functions
### gridSearch
This function performs a grid search to find the best performing combination of parameters (included are learning rate, number of observations in node, and maximum depth). The result from this function is a list object containing 1) as data.frame of the parameter combinations, the best performing number of trees, and their evaluation statistics and 2) a data.frame containing all of the information needed to construct a multipanel plot of deviance curves. 

### bootstrapGBM
This function evaluates a gbm model with a certain combination of parameters. The number of bootstrap iterations and whether the bootstrap run uses observed data or reshuffled labels on your data.frame ("null" bootstrap for evaluation statistic correction purposes) can be specified. The result of this function is a list object containing 1) a data.frame of the evaluation statistics pertaining to each bootstrap run (easy to use apply and get means from this element), 2) output used to create partial dependency plots using the partialPlot function described below, and 3) predictions on the dataset you entered. 

### partialPlot
This function creates partial dependency plots similar to Han et al. 2015 (Rodent reservoirs of future diseases). The result from this function is a multipanel plot showing the marginal effect on prediction on the left y-axis (and as the line plot layer) and frequency of the distribution of a trait on the right y-axis (as depicted by the histogram layer and understood by the breaks along the x-axis). It uses the output from an observed data run of the bootstrapGBM function and can be modified to only include variables of interest based on what is given to the vars argument. 
