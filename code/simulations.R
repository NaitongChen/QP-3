# load helper functions
source("helpers.R")

p = 20
n = 100
beta = rep(0, p)
beta[1] = 1
beta[2] = -1
beta[19] = 1
beta[20] = -1
signal = 1
mu_x = 0
sigma_x = 1
sigma_y = 1
a = 0.5

dat = generate_data(n, p, beta, signal, mu_x, sigma_x, sigma_y)
X = dat[[1]]
y = dat[[2]]