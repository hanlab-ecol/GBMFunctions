# GBMFunctions
A collection of functions designed to help in the parameterization, evaluation, and exploration of models using {gbm}.

## Current functions
### gridSearch
This function performs a grid search to find the best performing combination of parameters (included are learning rate, number of observations in node, and maximum depth). The result from this function is a list object containing 1) as data.frame of the parameter combinations, the best performing number of trees, and their evaluation statistics and 2) a data.frame containing all of the information needed to construct a multipanel plot of deviance curves. 

### spatialGridSearch
This function performs the same set of tasks as gridSearch, except that it allows for the the specification of groups to use for gbm's internal cross validation process. This requires that the gbmCrossVal function is sourced (this will add a line that allows for the specified groups to be used instead) and that a vector of group designations be supplied.

### bootstrapGBM
This function evaluates a gbm model with a certain combination of parameters. The number of bootstrap iterations and whether the bootstrap run uses observed data or reshuffled labels on your data.frame ("null" bootstrap for evaluation statistic correction purposes) can be specified. The result of this function is a list object containing 1) a data.frame of the evaluation statistics pertaining to each bootstrap run (easy to use apply and get means from this element), 2) output used to create partial dependency plots using the partialPlot function described below, and 3) predictions on the dataset you entered. 

### spatialBootstrap
This function, like spatialGridSearch, is the spatially explicit version of bootstrapGBM. This allows for the specification of specific groups to be used in gbm's internal cross validation process. As above, make sure that the gbmCrossVal function is loaded as you will likely run into an error without it (likely an object not found one, but who can predict these things sometimes?).

### partialPlot
This function creates partial dependency plots similar to Han et al. 2015 (Rodent reservoirs of future diseases). The result from this function is a multipanel plot showing the marginal effect on prediction on the left y-axis (and as the line plot layer) and frequency of the distribution of a trait on the right y-axis (as depicted by the histogram layer and understood by the breaks along the x-axis). It uses the output from an observed data run of the bootstrapGBM function and can be modified to only include variables of interest based on what is given to the vars argument. 

### group.block
This is a modified version of the group.blocks function from {ENMeval} to allow for groups larger than 4 to be specified. It will start assigning groups to the point furthest from the others. Do note that it is written to provide an equal number of background and presence points in each group, so the groupings may melt into one another (thanks to the background points). Will probably turn this into an argument to allow for discrete groupings. Requires {sp}.
