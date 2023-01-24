library(glmnet)
library(ggplot2)

generate_data <- function(n, p, beta, signal, mu_x, sigma_x, sigma_y) {
  X = matrix(rnorm(n*p, mu_x, sigma_x), nrow=n)
  y = X %*% (signal * beta) + rnorm(n, 0, sigma_y)
  
  return(list(X,y))
}

split_data <- function(X, y, a) {
  n = nrow(X)
  split_size = floor(n * a)
  shuffle_inds = sample(1:n, n)

  first_half <- shuffle_inds[1:split_size]
  second_half <- shuffle_inds[(split_size+1):n]
  
  return(list(X[first_half,], y[first_half], X[second_half,], y[second_half]))
}

fission_data_1 <- function(X, y, tau, sigma_y) {
  n = nrow(X)
  Z = rnorm(n, 0, sigma_y)
  
  first_half <- y + (tau * Z)
  second_half <- y - (1/tau * Z)
  
  return(list(X, first_half, X, second_half, Z))
}

fission_data_2 <- function(X, y, tau, sigma_y) {
  n = nrow(X)
  Z <- y + rnorm(n, 0, sqrt(tau) * sigma_y)
  
  first_half <- Z
  
  return(list(X, first_half, X, y, Z))
}

variable_selection <- function(X, y) {
  n = nrow(X)
  p = ncol(X)
  cvfit <- cv.glmnet(X, drop(y), type.measure = "mse", family = "gaussian", 
                     alpha = 1, nfolds = n, grouped = FALSE, intercept = FALSE)
  final_fit <- glmnet(X, drop(y), type.measure = "mse", family = "gaussian", 
                      alpha = 1, lambda = cvfit$lambda.1se, grouped = FALSE, 
                      intercept = FALSE)
  selected <- (1:p)[drop(final_fit$beta) != 0]

  return(selected)
}

build_CI_split <- function(X, y, sigma_y, selected, sig) {
  X = matrix(drop(X[, selected]), ncol = length(selected), byrow = FALSE)
  p = ncol(X)
  n = nrow(X)
  df = data.frame(X,y)

  beta_hat = lm(y~.-1, data = df)$coefficients
  
  normal_mat = tryCatch({
    f = solve(t(X) %*% X)
  }, error = function(err) {
    f = solve(t(X) %*% X + (1e-8) * diag(p))
    return(f)
  })
  
  cov = (normal_mat %*% t(X) %*% (sigma_y^2 * diag(n)) 
         %*% X %*% normal_mat)

  z = qnorm(1-(sig/2/length(selected)))

  CIs = matrix(rep(0, p*2), nrow = p)

  for (i in 1:p) {
    CIs[i,1] = beta_hat[i] - z * cov[i,i]
    CIs[i,2] = beta_hat[i] + z * cov[i,i]
  }

  return(list(beta_hat, CIs))
}

build_CI_1 <- function(X, y, tau, sigma_y, selected, sig) {
  X = matrix(drop(X[, selected]), ncol = length(selected), byrow = FALSE)
  p = ncol(X)
  n = nrow(X)
  df = data.frame(X,y)
  
  beta_hat = lm(y~.-1, data = df)$coefficients
  
  normal_mat = tryCatch({
    f = solve(t(X) %*% X)
  }, error = function(err) {
    f = solve(t(X) %*% X + (1e-8) * diag(p))
    return(f)
  })
  
  cov = (1 + tau^-2) * (normal_mat %*% t(X) %*% (sigma_y^2 * diag(n)) 
                        %*% X %*% normal_mat)
  
  z = qnorm(1-(sig/2/length(selected)))
  
  CIs = matrix(rep(0, p*2), nrow = p)
  
  for (i in 1:p) {
    CIs[i,1] = beta_hat[i] - z * cov[i,i]
    CIs[i,2] = beta_hat[i] + z * cov[i,i]
  }
  
  return(list(beta_hat, CIs))
}

build_CI_2 <- function(X, y, tau, sigma_y, selected, sig) {
  X = matrix(drop(X[, selected]), ncol = length(selected), byrow = FALSE)
  p = ncol(X)
  n = nrow(X)
  df = data.frame(X,y)
  
  beta_hat = lm(y~.-1, data = df)$coefficients
  
  normal_mat = tryCatch({
    f = solve(t(X) %*% X)
  }, error = function(err) {
    f = solve(t(X) %*% X + (1e-8) * diag(p))
    return(f)
  })
  
  cov = (tau/(tau + 1)) * (normal_mat %*% t(X) %*% (sigma_y^2 * diag(n)) 
                           %*% X %*% normal_mat)
  
  z = qnorm(1-(sig/2/length(selected)))
  
  CIs = matrix(rep(0, p*2), nrow = p)
  
  for (i in 1:p) {
    CIs[i,1] = beta_hat[i] - z * cov[i,i]
    CIs[i,2] = beta_hat[i] + z * cov[i,i]
  }
  
  return(list(beta_hat, CIs))
}

plot_single_trial <- function(trial, method, n, sigma_y, beta, signal, 
                              selected, beta_hat, CIs) {
  p = length(beta)
  index <- 1:p
  lowers = rep(0, p)
  uppers = rep(0, p)
  beta_hats = rep(0, p)
  lowers[selected] = CIs[,1]
  uppers[selected] = CIs[,2]
  beta_hats[selected] = beta_hat
  beta = signal * beta
  
  df = data.frame(index, beta, beta_hats, lowers, uppers)
  df_selected = df[selected,]
  
  colors = c()
  for (i in 1:length(selected)) {
    if (CIs[i,1] <= df_selected$beta[i] & CIs[i,2] >= df_selected$beta[i]) {
      colors = c(colors, "blue") 
    } else {
      colors = c(colors, "red") 
    }
  }
  
  ggplot(df_selected, aes(index, beta)) +
    geom_point(shape=4, size=4, color="blue") +
    geom_errorbar(aes(ymin = lowers, ymax = uppers), color=colors) +
    geom_point(aes(index, beta_hats), color=colors) +
    geom_point(aes(index, beta), df)

  ggsave(paste0("../plot/", method, "_", as.character(n), "_", 
                as.character(sigma_y), "_", as.character(trial), ".png"))
}

metric_single_trial <- function(beta, signal, selected, beta_hat, CIs) {
  p = length(beta)
  index <- 1:p
  lowers = rep(0, p)
  uppers = rep(0, p)
  beta_hats = rep(0, p)
  lowers[selected] = CIs[,1]
  uppers[selected] = CIs[,2]
  beta_hats[selected] = beta_hat
  beta = signal * beta
  
  df = data.frame(index, beta, beta_hats, lowers, uppers)
  df_selected = df[selected,]
  
  pow = length(intersect((1:p)[beta != 0], selected)) / length((1:p)[beta != 0])
  prec = length(intersect((1:p)[beta != 0], selected)) / length(selected)
  avg_CI_length = mean(CIs[,2] - CIs[,1])
  
  FCR = 0
  for (i in 1:length(selected)) {
    if (CIs[i,1] <= df_selected$beta[i] & CIs[i,2] >= df_selected$beta[i]) {
      FCR = FCR + 0
    } else {
      FCR = FCR + 1
    }
  }
  FCR = FCR / max(1, length(selected))
  
  bhat_full = rep(0, p)
  bhat_full[selected] = beta_hat
  L2_err = norm(as.matrix(signal*beta - bhat_full), type="2")
  
  return(list(pow, prec, avg_CI_length, FCR, L2_err))
}

plot_metric <- function(trial, ns, sigma_y, metric, dat1, dat2, dat3) {
  n_size = length(ns)
  size = ns
  medians = rep(0, n_size)
  first_q = rep(0, n_size)
  third_q = rep(0, n_size)
  for (i in 1:n_size) {
    medians[i] = median(dat1[[i]])
    first_q[i] = quantile(dat1[[i]], 0.25)
    third_q[i] = quantile(dat1[[i]], 0.75)
  }
  
  df1 = data.frame(size, medians, first_q, third_q)
  
  medians = rep(0, n_size)
  first_q = rep(0, n_size)
  third_q = rep(0, n_size)
  for (i in 1:n_size) {
    medians[i] = median(dat2[[i]])
    first_q[i] = quantile(dat2[[i]], 0.25)
    third_q[i] = quantile(dat2[[i]], 0.75)
  }
  
  df2 = data.frame(size, medians, first_q, third_q)
  
  medians = rep(0, n_size)
  first_q = rep(0, n_size)
  third_q = rep(0, n_size)
  for (i in 1:n_size) {
    medians[i] = median(dat3[[i]])
    first_q[i] = quantile(dat3[[i]], 0.25)
    third_q[i] = quantile(dat3[[i]], 0.75)
  }
  
  df3 = data.frame(size, medians, first_q, third_q)
  
  ggplot() +
    geom_point(aes(size, medians), df1, color="red") +
    #geom_errorbar(aes(x = df1$size, ymin = df1$first_q, ymax = df1$third_q), color="red") +
    geom_point(aes(size, medians), df2, color="blue") +
    #geom_errorbar(aes(x = df2$size, ymin = df2$first_q, ymax = df2$third_q), color="blue") +
    geom_point(aes(size, medians), df3, color="green") +
    #geom_errorbar(aes(x = df3$size, ymin = df3$first_q, ymax = df3$third_q), color="green") +
    ylab(metric)
  
  ggsave(paste0("../plot/", metric, "_", 
                as.character(sigma_y), ".png"))
}