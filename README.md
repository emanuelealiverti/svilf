# Svilf: Stratified stochastic variational inference for high-dimensional network factor model	

This repository is associated with the article: Aliverti, E. and Russo, M. (2020). Stratified stochastic variational inference for high-dimensional network factor model	[[link](https://arxiv.org/abs/0000.000)].

> In this article, we focus on the latent factor model for networks, a very general
> approach which encompasses the latent distance model and the stochastic block-model as special cases. We
> develop scalable algorithms to conduct approximate Bayesian inference via stochastic optimization. Leveraging
> sparse representations of network data, the proposed algorithms show massive computational and storage
> benefits, and allow to conduct inference in settings with thousands of nodes.



The repository contains the `R` package implementing `svilf`.
We recommend installing `svilf` locally in the `lib/` folder, for example with

```
R CMD install -l lib ./svilf_0.1.tar.gz
```

The package provides wrappers around core `c++` implementations, relying on Armdadillo's `SpMat` [format](https://cran.r-project.org/web/packages/RcppArmadillo/vignettes/RcppArmadillo-sparseMatrix.pdf).
Wrappers are written in `R` and provide easy tools to specify different settings. We also provide default parameter implementation, which in our experience works well in many settings.

The scripts in the folders `APPLICATION` and `SIMULATIONS` reproduce results from the paper, and are are set-up to use `./lib` as the local library folder.
