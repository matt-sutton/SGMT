# Temporary package for implementing Sparse Group Multi-Task methods for pleiotropy 

This is a temporary package for the implementation of penalised multi-task methods that incorporate sparsity in effect estimation for multiple logistic regression studies. The package builds on some of the structured R code from glmnet [Noah, 2011] for cross validation (with reference and citation in the code headers see cv_method(s) in the R directory) and extends this to allow for analysis of multiple studies. The underlying algorithm is implemented using C++ and RcppNumerical for optimisation of the logistic likelihood. 

Variables are selected across studies at an individual, group or sparse group structure. This is facilitated through the use of penalisation methods and these methods are implemented using the alternating direction method of multipliers (ADMM) algorithm [Boyd 2011]. 

The current release allows for:

1. Implementation of variable selection via sparse, group or group sparse penalisation for multiple studies with logistic regression likelihoods. If a variable or group of variables are selected, they are selected across all studies.

2. Cross validation methods for choosing the tuning parameters in all penalisation methods.

3. Example code that includes illustrative bootstrapping code to detect variability in selected variables/groups of variables. 

4. Choosing adaptive weights for the penalisation in the adaptive versions of the method.

Future versions may build out on the flexibility of likelihoods offered across studies or alternative penalisation methods. 

## Installation:

You can install the package using devtools:
```{r,}
library(devtools)
install_github("matt-sutton/SGMT")
```

## Example Code:

The plot below was generated from the example_simulation_run.R code which replicates a single run of the methods on a dataset from scenario 2 of the simulation study from the paper. The top row is the true regression coefficients used to simulate the datasets, other rows correspond to different methods.

![Scenario2](https://github.com/matt-sutton/SGMT/blob/main/Scenario2.png)

These results are the estimated coefficients for the selected variables (using bootstrapping as described in the paper) and FDR control for ASSET. 


## References:

Simon, Noah, Jerome Friedman, Trevor Hastie, and Rob Tibshirani. 2011.
“Regularization Paths for Cox’s Proportional Hazards Model via
Coordinate Descent.” *Journal of Statistical Software, Articles* 39 (5): 1–13. <https://doi.org/10.18637/jss.v039.i05>.

Stephen Boyd, Neal Parikh, Eric Chu., Borja Peleato., Jonathan Eckstein. 2011.
"Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers."
*Foundations and Trendsin Machine Learning* 3(1), 1–122.
  
