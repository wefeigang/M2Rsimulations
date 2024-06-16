# code to implement cross-validation for graphical lasso

library(glasso)
library(MASS) # For mvrnorm function
library(CVglasso) # to test

set.seed(13)

cv_glasso <- function(X, nlam = 100, K=10, tol = 1e-04, maxit = 10000){
  S <- cov(X)
  diag(S) <- 0
  lam.max <- max(abs(S))
  lam.min <- ifelse(nrow(X) > ncol(X), 10^-4, 0.01)
  lambda_seq <- 10^seq(log10(lam.min), log10(lam.max), length = nlam)
  # Number of observations
  n <- nrow(X)
  
  # Shuffle indices to ensure randomised folds
  ind <- sample(n)
  
  # Initialise a matrix to store cross-validation errors
  cv_errors <- matrix(0, nrow = length(lambda_seq), ncol = K)
  
  # Perform K-fold cross-validation
  for (k in 1:K) {
    # Determine the indices for the validation set
    val_ind <- ind[(1 + floor((k - 1) * n/K)):floor(k * n/K)]
    
    # Separate training and validation sets
    train_data <- X[-val_ind, , drop = FALSE] # These include whole rows
    val_data <- X[val_ind, , drop = FALSE]
    
    # Centre the training and validation data (just to "scale"/centralise the data so do I need to scale data later when computing metrics?)
    X_bar <- colMeans(train_data)
    train_data <- scale(train_data, center = X_bar, scale = FALSE)
    val_data <- scale(val_data, center = X_bar, scale = FALSE)
    
    # Sample covariance matrices
    S_train <- cov(train_data)
    S_val <- cov(val_data)
    
    # Loop over all lambda values
    for (i in seq_along(lambda_seq)){
      lambda <- lambda_seq[i]
      
      # Fit graphical lasso on training data (what is penalize diagonal?)
      glasso_fit <- glasso(S_train, rho = lambda, thr = tol, maxit = maxit, penalize.diagonal = FALSE)
      
      # Estimated precision matrix
      est_prec <- glasso_fit$wi
      
      # Compute the observed negative validation log-likelihood (need citation found from github)
      log_det_term <- determinant(est_prec, logarithm = TRUE)$modulus
      trace_term <- sum(diag(est_prec %*% S_val))
      cv_errors[i, k] <- 0.5 * nrow(X) * (trace_term - log_det_term) # is nrow(X) necessary?
    }
  }
  
  # Average CV error for each lambda by finding mean across all k (row)
  mean_cv_errors <- rowMeans(cv_errors)
  
  # Find the lambda with the minimum CV error
  best_lambda <- lambda_seq[which.min(mean_cv_errors)]
  
  return(list(best_lambda = best_lambda, cv_errors = cv_errors))
}

# example

# Make sparse precision matrix
sparse_matrix <- function(p){
  true_prec <- diag(1, p)
  #true_prec[abs(row(true_prec) - col(true_prec)) == 1] <- 0.5
  true_prec[abs(row(true_prec) - col(true_prec)) == 1] <- 0.4
  true_prec[abs(row(true_prec) - col(true_prec)) == 2] <- 0.2
  true_prec[abs(row(true_prec) - col(true_prec)) == 3] <- 0.2
  return(true_prec)
}

set.seed(123)
p <- 20
n <- 10
true_prec <- sparse_matrix(p)
true_cov <- solve(true_prec)
data <- mvrnorm(n = n, numeric(p), Sigma = true_cov)
data <- scale(data)

result <- cv_glasso(X = data, K=5)

print(result$cv_errors)
print(result$best_lambda)

cv_glasso_result <- CVglasso(X = data)
print(cv_glasso_result$Tuning)
#print(cv_glasso_result$CV.error)
