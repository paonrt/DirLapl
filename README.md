# DirLapl package
R package for posterior computation under a _Dirichlet-Laplace prior_.

For the linear regression model and the logistic regression model, customized
Gibbs samplers are available. For a generic likelihood, it is used an elliptical
slice sampling.

Required dependencies:

* BayesLogit
* GIGrvg
* glmnet
* statmod

Install the package with:

```R
devtools::install_github("paonrt/DirLapl")
```
