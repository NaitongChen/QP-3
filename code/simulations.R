# load helper functions
setwd("C:/projects/Qualifying papers/QP-3/code")
source("helpers.R")
library(torch)

# set parameters
p = 20
ns = c(10, 20, 50, 100)
beta = rep(0, p)
beta[1] = 1
beta[2] = -1
beta[19] = 1
beta[20] = -1
signal = 1
mu_x = 0
sigma_x = 1
sigma_ys = c(0.25, 0.5, 1, 2)
a = 0.5
sig = 0.05
n_trial = 20
methods = c("split", "p1", "p2")
metrics = c("power", "precision", "len", "FCR", "L2")
tau = 1

# storing output
pow_split = torch_zeros(c(length(ns), length(sigma_ys), n_trial)) # n by sigma_y by n_trial
prec_split = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
len_split = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
fcr_split = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
l2_split = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
noselect_split = torch_zeros(c(length(ns), length(sigma_ys), n_trial))

pow_p1 = torch_zeros(c(length(ns), length(sigma_ys), n_trial)) # n by sigma_y
prec_p1 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
len_p1 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
fcr_p1 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
l2_p1 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
noselect_p1 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))

pow_p2 = torch_zeros(c(length(ns), length(sigma_ys), n_trial)) # n by sigma_y
prec_p2 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
len_p2 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
fcr_p2 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
l2_p2 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))
noselect_p2 = torch_zeros(c(length(ns), length(sigma_ys), n_trial))

# run simulations
for (i in 1:n_trial) {
  print(paste0(as.character(i), "/", as.character(n_trial)))
  set.seed(i)
  for (j in 1:length(ns)) {
    print(paste0("n=", as.character(ns[j])))
    for (k in 1:length(sigma_ys)) {
      print(paste0("sigma_y=", as.character(sigma_ys[k])))
      dat = generate_data(ns[j], p, beta, signal, mu_x, sigma_x, sigma_ys[k])
      X = dat[[1]]
      y = dat[[2]]
      true_mu = dat[[3]]
      
      # split
      dat = split_data(X, y, a)
      X1 = dat[[1]]
      y1 = dat[[2]]
      X2 = dat[[3]] 
      y2 = dat[[4]]
      idx = dat[[5]]
      selected = variable_selection(X1, y1)
      if (length(selected) == 0) {
        noselect_split[j,k,i] = 1
      } else {
        dat = build_CI_split(X2, y2, sigma_ys[k], selected, sig)
        beta_hat = dat[[1]]
        CIs = dat[[2]]
        plot_single_trial(i, methods[1], ns[j], sigma_ys[k], beta, signal, 
                          selected, beta_hat, CIs)
        dat = metric_single_trial(beta, signal, selected, 
                                  beta_hat, CIs, X2, y2, true_mu[idx])
        pow_split[j,k,i] = dat[[1]]
        prec_split[j,k,i] = dat[[2]]
        len_split[j,k,i] = dat[[3]]
        fcr_split[j,k,i] = dat[[4]]
        l2_split[j,k,i] = dat[[5]]
      }
      
      # p1
      dat = fission_data_1(X, y, tau, sigma_ys[k])
      X1 = dat[[1]]
      y1 = dat[[2]]
      X2 = dat[[3]] 
      y2 = dat[[4]]
      selected = variable_selection(X1, y1)
      if (length(selected) == 0) {
        noselect_p1[j,k,i] = 1
      } else {
        dat = build_CI_1(X2, y2, tau, sigma_ys[k], selected, sig)
        beta_hat = dat[[1]]
        CIs = dat[[2]]
        plot_single_trial(i, methods[2], ns[j], sigma_ys[k], beta, signal, 
                          selected, beta_hat, CIs)
        dat = metric_single_trial(beta, signal, selected, 
                                  beta_hat, CIs, X2, y2, true_mu)
        pow_p1[j,k,i] = dat[[1]]
        prec_p1[j,k,i] = dat[[2]]
        len_p1[j,k,i] = dat[[3]]
        fcr_p1[j,k,i] = dat[[4]]
        l2_p1[j,k,i] = dat[[5]] 
      }
      
      # p2
      dat = fission_data_2(X, y, tau, sigma_ys[k])
      X1 = dat[[1]]
      y1 = dat[[2]]
      X2 = dat[[3]] 
      y2 = dat[[4]]
      selected = variable_selection(X1, y1)
      if (length(selected) == 0) {
        noselect_p2[j,k,i] = 1
      } else {
        dat = build_CI_2(X2, y2, tau, sigma_ys[k], selected, sig)
        beta_hat = dat[[1]]
        CIs = dat[[2]]
        plot_single_trial(i, methods[3], ns[j], sigma_ys[k], beta, signal, 
                          selected, beta_hat, CIs)
        dat = metric_single_trial(beta, signal, selected, 
                                  beta_hat, CIs, X2, y2, true_mu)
        pow_p2[j,k,i] = dat[[1]]
        prec_p2[j,k,i] = dat[[2]]
        len_p2[j,k,i] = dat[[3]]
        fcr_p2[j,k,i] = dat[[4]]
        l2_p2[j,k,i] = dat[[5]] 
      }
    }
  }
}

save.image(file="end_of_simulations.RData")

for (k in 1:length(sigma_ys)) {
  # split
  dat1 = list()
  for (i in 1:length(ns)) {
    dat1[[i]] =  as_array(pow_split[i,k,])
  }
  
  dat2 = list()
  for (i in 1:length(ns)) {
    dat2[[i]] = as_array(pow_p1[i,k,])
  }
  
  dat3 = list()
  for (i in 1:length(ns)) {
    dat3[[i]] = as_array(pow_p2[i,k,])
  }
  
  plot_metric(trial, ns, sigma_ys[k], metrics[1], dat1, dat2, dat3)
  
  dat1 = list()
  for (i in 1:length(ns)) {
    dat1[[i]] = as_array(prec_split[i,k,])
  }
  
  dat2 = list()
  for (i in 1:length(ns)) {
    dat2[[i]] = as_array(prec_p1[i,k,])
  }
  
  dat3 = list()
  for (i in 1:length(ns)) {
    dat3[[i]] = as_array(prec_p2[i,k,])
  }
  
  plot_metric(trial, ns, sigma_ys[k], metrics[2], dat1, dat2, dat3)
  
  dat1 = list()
  for (i in 1:length(ns)) {
    idx = (1:n_trial)[as_array(noselect_split[i,k,]) == 0]
    dat1[[i]] = as_array(len_split[i,k,][idx])
  }
  
  dat2 = list()
  for (i in 1:length(ns)) {
    idx = (1:n_trial)[as_array(noselect_p1[i,k,]) == 0]
    dat2[[i]] = as_array(len_p1[i,k,][idx])
  }
  
  dat3 = list()
  for (i in 1:length(ns)) {
    idx = (1:n_trial)[as_array(noselect_p2[i,k,]) == 0]
    dat3[[i]] = as_array(len_p2[i,k,][idx])
  }
  
  plot_metric(trial, ns, sigma_ys[k], metrics[3], dat1, dat2, dat3)
  
  dat1 = list()
  for (i in 1:length(ns)) {
    dat1[[i]] = as_array(fcr_split[i,k,])
  }
  
  dat2 = list()
  for (i in 1:length(ns)) {
    dat2[[i]] = as_array(fcr_p1[i,k,])
  }
  
  dat3 = list()
  for (i in 1:length(ns)) {
    dat3[[i]] = as_array(fcr_p2[i,k,])
  }
  
  plot_metric(trial, ns, sigma_ys[k], metrics[4], dat1, dat2, dat3)
  
  dat1 = list()
  for (i in 1:length(ns)) {
    idx = (1:n_trial)[as_array(noselect_split[i,k,]) == 0]
    dat1[[i]] = as_array(l2_split[i,k,][idx])
  }
  
  dat2 = list()
  for (i in 1:length(ns)) {
    idx = (1:n_trial)[as_array(noselect_p1[i,k,]) == 0]
    dat2[[i]] = as_array(l2_p1[i,k,][idx])
  }
  
  dat3 = list()
  for (i in 1:length(ns)) {
    idx = (1:n_trial)[as_array(noselect_p2[i,k,]) == 0]
    dat3[[i]] = as_array(l2_p2[i,k,][idx])
  }
  
  plot_metric(trial, ns, sigma_ys[k], metrics[5], dat1, dat2, dat3)
}