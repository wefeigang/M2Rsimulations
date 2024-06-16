library(glasso) # glasso package
library(CVglasso) # extension of glasso package which uses cross-validation
library(MASS) # contains mvrnorm
library(progress) # for progress bar

# Make sparse precision matrix
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

# Cross validation
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

# Simulation function
glasso_simulation <- function(n, p, k, data){
  # Forming a tridiagonal precision matrix
  #true_prec <- sparse_matrix(p = p, k = k)
  #true_cov <- solve(true_prec) # round to 5 d.p. to get rid of slight inaccuracies when inverting matrix
  #data <- mvrnorm(n = n, numeric(p), Sigma = true_cov)
  # data <- scale(data)
  # should I scale the data?
  cv_lambda <- cv_glasso(X = data)$best_lambda
  est_prec <- glasso(cov(data), rho = cv_lambda)$wi
  est_cov <- glasso(cov(data), rho = cv_lambda)$w
  return(list(est_cov = est_cov, est_prec = est_prec))
}

# Kullback–Leibler divergence for matrices note there is no mean terms as they are equal
kl_divergence <- function(est_matrix, true_matrix){
  p <- ncol(true_matrix)
  log_det_ratio <- log(det(est_matrix)) - log(det(true_matrix))
  trace_term <- sum(diag(solve(est_matrix) %*% true_matrix))
  return(0.5 * (trace_term - p + log_det_ratio))
}

# Finding number of correct non-zero entries/total non-zero entries of estimated precision matrix
metric1 <- function(est_matrix, true_matrix){
  num_correct <- sum(est_matrix != 0 & true_matrix != 0)
  total <- sum(est_matrix != 0)
  
  return(num_correct/total)
}

# Finding number of correct zero entries/total zero entries of estimated precision matrix
metric2 <- function(est_matrix, true_matrix){
  num_correct <- sum(est_matrix == 0 & true_matrix == 0)
  total <- sum(est_matrix == 0)
  
  return(num_correct/total)
}

# To find the maximum difference in the entries using frobenius norm/ number of entries (on precision)
evaluate_frob <- function(est_matrix, true_matrix, p){
  diff_mat <- est_matrix - true_matrix
  return(norm(diff_mat, type = "F")/norm(true_matrix, type = "F"))
}

# F-score (useful when sparse)
evaluate_F <- function(estimated_precision, true_precision){
  tp <- sum(true_precision != 0 & estimated_precision != 0)
  tn <- sum(true_precision == 0 & estimated_precision == 0)
  fp <- sum(true_precision == 0 & estimated_precision != 0)
  fn <- sum(true_precision != 0 & estimated_precision == 0)
  
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)
  f1_score <- 2 * (precision * recall) / (precision + recall)
  
  return(f1_score)
}

# Define the range of values for n and p
n_values <- c(20, 40, 60, 100)
p_values <- c(10, 20, 30, 50)
#n_values <- c(20, 30)
#p_values <- c(10, 20)

n_iter <- 3 * length(n_values) * length(p_values)

# Progress bar
pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                       total = n_iter,
                       complete = "=",   # Completion bar character
                       incomplete = "-", # Incomplete bar character
                       current = ">",    # Current bar character
                       clear = FALSE,    # If TRUE, clears the bar when finish
                       width = 100)      # Width of the progress bar
              
# Initialize an empty data frame to store the results
results <- data.frame(n = integer(), p = integer(), mean_kl_divergence = numeric())

start.time <- Sys.time() # Start time to measure execution time

# Loop through each combination of n and p and for each type of matrix
for (k in 1:3){
  for (n in n_values){
    for (p in p_values){
      pb$tick() # progress bar tick
      
      # Generate the true covariance matrix
      true_prec <- sparse_matrix(p = p, k = k)
      true_cov <- solve(true_prec)
      data <- mvrnorm(n = n, numeric(p), Sigma = true_cov)
      
      # Replicate the simulation multiple times and compute the mean KL divergence
      kl_scores <- replicate(10, kl_divergence(glasso_simulation(n = n, p = p, k = k, data = data)$est_cov, true_matrix = true_cov))
      kl_mean <- mean(kl_scores)
      #metric_values1 <- replicate(10, metric1(glasso_simulation(n = n, p = p, k = k), true_prec))
      #metric_values2 <- replicate(10, metric2(glasso_simulation(n = n, p = p, k = k), true_prec))
      #mean1 <- mean(metric_values1)
      #mean2 <- mean(metric_values2)
      #mean <- (mean1 + mean2)/2
      
      frob_scores <- replicate(10, evaluate_frob(glasso_simulation(n = n, p = p, k = k, data = data)$est_prec, true_prec, p = p))
      frob_mean <- mean(frob_scores)
      
      f_scores <- replicate(10, evaluate_F(estimated_precision = glasso_simulation(n = n, p = p, k = k, data = data)$est_prec, true_precision = true_prec))
      f_mean <- mean(f_scores)
      
      # Store the results in the data frame
      results <- rbind(results, data.frame(n = n, p = p, k = k, f_score = f_mean, KL_divergence = kl_mean, frobenius = frob_mean))
    }
  }
}

end.time <- Sys.time() # End time
time.taken <- round(end.time - start.time,2)

# Print the results as a table
print(results)
print(time.taken)
