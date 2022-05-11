# Installing the library

Install the package provided in this folder following the instruction on
`README.md` or running

``` r
install.packages("svilf_0.1.tar.gz",lib = "lib", type = "source", repos = NULL)
```

``` r
library(svilf, lib = "lib")
args(svilf)
```

    ## function (e, H = 2, intercept = T, eigen_init = F, prop = 2, 
    ##     link = "logit", get_samples = 0, opts = svilf_options()) 
    ## NULL

The main function `svilf` takes an edge-list as input.

We generate an artificial network via `igraph`, and convert it into ad
edge-list. This is the format which is naturally used by `svilf`

``` r
library(igraph)
set.seed(1)
net = sample_dot_product(sample_sphere_surface(dim = 2, n = 100))
el = get.edgelist(net)
```

``` r
# defualt options
m_log = svilf(el)

# probit link, setting H = 3
m_prob = svilf(el, H = 3, link = "probit")

# Adaptive sampling (denoted as ADA in the paper)
m_logADA = svilf(el, opts = svilf_options(sample_adaptive = T))
m_probADA = svilf(el, link = "probit", opts = svilf_options(sample_adaptive = T))

# reduce tolerance
m_log = svilf(el, opts = svilf_options(tol = 1e-4))

# sensible initialization setting Z to eigenvectors
m_log = svilf(el, eigen_init=T)
```
