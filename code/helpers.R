library(glmnet)

generate_data <- function(n, p, beta, signal, mu_x, sigma_x, sigma_y) {
  X = matrix(rnorm(n*p, mu_x, sigma_x), nrow=n)
  y = X %*% (signal * beta) + rnorm(n, 0, sigma_y)
  
  return(list(X,y))
}

split_data <- function(X, y, a) {
  n = nrow(X)
  split_size = floor(n * a)
  shuffle_inds = sample(1:n, n)

  first_half = shuffle_inds[1:split_size]
  second_half = shuffle_inds[split_size+1:n]
  
  return(list(X[first_half,], y[fist_half], X[second_half,], y[second_half]))
}

fission_data_1 <- function(X, y, tau, sigma_y) {
  n = nrow(X)
  Z = rnorm(n, 0, sigma_y)
  
  first_half = y + (tau * Z)
  second_half = y - (1/tau * Z)
  
  return(list(X, fist_half, X, second_half, Z))
}

fission_data_2 <- function(X, y, tau, sigma_y) {
  n = nrow(X)
  Z = y + rnorm(n, 0, sqrt(tau) * sigma_y)
  
  first_half = Z
  
  return(list(X, fist_half, X, y, Z))
}

variable_selection <- function(X, y) {
  n = nrow(X)
  p = ncol(X)
  cvfit <- cv.glmnet(X, drop(y), type.measure = "mse", family = "gaussian", 
                     alpha = 1, nfolds = n, grouped = FALSE)
  final_fit <- glmnet(X, drop(y), type.measure = "mse", family = "gaussian", 
                      alpha = 1, lambda = cvfit$lambda.1se)
  selected <- (1:p)[drop(final_fit$beta) != 0]
  
  return(selected)
}

build_CI_split <- function(X, y, sigma_y, selected) {
  # todo
}

build_CI_1 <- function(X, y, tau, sigma_y, selected) {
  # todo
}

build_CI_2 <- function(X, y, tau, sigma_y, selected) {
  # todo
}