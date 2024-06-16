library(glasso)
library(MASS)
library(CVglasso)


sparse_matrix <- function(p, k){ # k is the number from paper by Luan AR(k)
  true_prec <- diag(1, p)
  if (k == 1){
    true_prec[abs(row(true_prec) - col(true_prec)) == 1] <- 0.5
  }
  if (k == 2){
    true_prec[abs(row(true_prec) - col(true_prec)) == 1] <- 0.5
    true_prec[abs(row(true_prec) - col(true_prec)) == 2] <- 0.25
  }
  if (k == 3){
    true_prec[abs(row(true_prec) - col(true_prec)) == 1] <- 0.4
    true_prec[abs(row(true_prec) - col(true_prec)) == 2] <- 0.2
    true_prec[abs(row(true_prec) - col(true_prec)) == 3] <- 0.2
  }
  return(true_prec)
}

# Define the cv_glasso function as provided
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

  return(list(best_lambda = best_lambda, cv_errors = cv_errors, mean_cv_errors = mean_cv_errors, lambda_seq = lambda_seq))
}


# Plot
par(mfrow=c(1, 3), mar=c(5,5,5,5))

for (k in 1:3){
  # Generate a sample dataset (you can replace this with your actual data)
  set.seed(123)
  n = 1000
  p = 20
  true_prec <- sparse_matrix(p, k = k)
  true_cov <- solve(true_prec) 
  X <- mvrnorm(n = n, numeric(p), Sigma = true_cov)
  
  # Calculate cross-validation errors
  results <- cv_glasso(X, nlam = 100)
  
  # Plotting the cross-validation errors against lambda values
  lambda_seq <- results$lambda_seq
  mean_cv_errors <- results$mean_cv_errors
  best_lambda <- results$best_lambda
  
  plot(lambda_seq, mean_cv_errors, main = paste0("Model ", k), pch = 20, type="o",xlab=expression(lambda), ylab="CV error", lwd = 2, cex = 2, cex.axis = 2, cex.lab = 2, cex.main = 2, log = "x")
  
  abline(v = best_lambda, col = "red", lty = 2, lwd=4)
  
  # Add custom text for the best lambda to the right of the vertical line
  text_x <- best_lambda  # Adjust this value as needed to move the text to the right
  text_y <- 1/2 * (max(mean_cv_errors) - min(mean_cv_errors)) + min(mean_cv_errors) # Place the text in the middle of the y-axis
  
  # Add text to indicate the best lambda at specified coordinates
  text(text_x, text_y, labels = bquote(lambda == .(round(best_lambda, 3))), col = "red", cex = 1.5, pos = 4)
}

