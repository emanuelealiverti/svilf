require(rstan)
tmp = stan_model("./factor.stan")
save.image("factor.RData")
