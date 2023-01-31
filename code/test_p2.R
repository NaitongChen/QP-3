# load helper functions
setwd("C:/projects/Qualifying papers/QP-3/code")
source("helpers.R")
library(torch)

# set parameters
p = 20
beta = rep(0, p)
beta[1] = 1
beta[2] = -1
beta[19] = 1
beta[20] = -1
signal = 1
mu_x = 0
sigma_x = 1
a = 0.5
sig = 0.05
n_trial = 100
methods = c("split", "p1", "p2")
metrics = c("power", "precision", "len", "FCR", "L2")
tau = 1

sigma_y = 1
n = 10

set.seed(runif(1, 1, 1000))
dat = generate_data(n, p, beta, signal, mu_x, sigma_x, sigma_y)
X = dat[[1]]
y = dat[[2]]
true_mu = dat[[3]]

# p1
dat = fission_data_2(X, y, tau, sigma_y)
X1 = dat[[1]]
y1 = dat[[2]]
X2 = dat[[3]] 
y2 = dat[[4]]
selected = variable_selection(X1, y1)
if (length(selected) > 0) {
  dat = build_CI_2(X2, y2, tau, sigma_y, selected, sig)
  beta_hat = dat[[1]]
  CIs = dat[[2]]
  
  plot_single_trial(1, methods[3], n, sigma_y, beta, signal,
                    selected, beta_hat, CIs, true_mu, X2)
  
  dat = metric_single_trial(beta, signal, selected, 
                            beta_hat, CIs, X2, y2, true_mu)
  pow_split = dat[[1]]
  prec_split = dat[[2]]
  len_split = dat[[3]]
  fcr_split = dat[[4]]
  l2_split = dat[[5]]
  
  print(pow_split)
  print(prec_split)
  print(len_split)
  print(fcr_split)
  print(l2_split)
}