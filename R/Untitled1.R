###############################################################################
# Load required package for multivariate normal simulation
###############################################################################
library(MASS)
library(Matrix)  # for bdiag if needed

###############################################################################
# 1) Data Generation Function
###############################################################################
simulate_data <- function(n1, m, q, q_a, r,
                          rhoX = 0, rhoM = 0,
                          sigma1 = 1,  # new: factor for mediator noise
                          sigma2 = 1,  # new: factor for outcome noise
                          pathway = c("complete","partial","none")) {
  
  pathway <- match.arg(pathway)
  
  ## 1) Generate exposures X (n1 x m) using compound-symmetry correlation with parameter rhoX
  SigmaX <- (1 - rhoX) * diag(m) + rhoX * matrix(1, nrow = m, ncol = m)
  X_unscaled <- mvrnorm(n1, mu = rep(0, m), Sigma = SigmaX)
  X <- scale(X_unscaled)  # standardize columns
  
  ## 2) Build block covariance for mediators M:
  #    - The first q_a mediators have correlation = rhoM.
  #    - The remaining (q - q_a) mediators are independent.
  q_inactive <- q - q_a
  block_active <- (1 - rhoM) * diag(q_a) + rhoM * matrix(1, nrow = q_a, ncol = q_a)
  if(q_inactive > 0){
    block_inactive <- diag(q_inactive)
    SigmaM <- as.matrix(bdiag(block_active, block_inactive))
  } else {
    SigmaM <- as.matrix(block_active)
  }
  
  # Multiply by sigma1^2 so that the noise has variance sigma1^2 (and correlations remain the same)
  SigmaM_scaled <- sigma1^2 * SigmaM
  
  # Generate mediator noise
  E_M <- mvrnorm(n1, mu = rep(0, q), Sigma = SigmaM_scaled)
  
  ## 3) Specify true effect parameters
  # Initialize loading matrix (for X -> M), direct effect vector gamma, and mediator effect vector eta
  alpha_mat <- matrix(0, nrow = m, ncol = q)
  gamma_vec <- rep(0, m)
  eta_vec   <- rep(0, q)
  
  # Active indices (for exposures and mediators)
  active_expos <- seq_len(r)
  active_mediators <- seq_len(q_a)
  
  if (pathway == "none") {
    # No mediation: all effect is direct.
    gamma_vec[active_expos] <- 1
    # alpha_mat and eta_vec remain all zeros.
    
  } else if (pathway == "complete") {
    # Purely indirect: no direct effect.
    gamma_vec[active_expos] <- 0
    for (i in seq_len(r)) {
      if (i <= q_a) {  # only assign for mediators with signal
        alpha_mat[i, i] <- 1
        eta_vec[i] <- 1
      }
    }
    
  } else if (pathway == "partial") {
    # Partial mediation: direct effect is 0.5 and indirect effect (via mediators) is 0.5.
    gamma_vec[active_expos] <- 1
    for (i in seq_len(r)) {
      if (i <= q_a) {
        alpha_mat[i, i] <- 0.1
        eta_vec[i] <- 0.1
      }
    }
  }
  
  ## 4) Generate mediators: M = X %*% alpha_mat + E_M
  M_unscaled <- X %*% alpha_mat + E_M
  M <- scale(M_unscaled)
  
  ## 5) Generate outcome: Y = X %*% gamma_vec + M %*% eta_vec + noise
  eps <- rnorm(n1, sd = sigma2)  # outcome noise with std dev = sigma2
  Y_raw <- as.vector(X %*% gamma_vec) + as.vector(M %*% eta_vec) + eps
  Y <- scale(Y_raw)
  Y <- as.vector(Y)
  
  list(X = X, M = M, Y = Y,
       active_expos = active_expos, active_mediators = active_mediators,
       true_alpha = alpha_mat, true_gamma = gamma_vec, true_eta = eta_vec)
}

###############################################################################
# 2) ADMM (same as before)
###############################################################################
optimize_weights <- function(X, M, Y, lambda_n, lambda_a, lambda_b, 
                             max_iter=50, tol=1e-4) {
  n <- nrow(X); m <- ncol(X); q <- ncol(M)
  
  # Initialize a and b (e.g., based on correlation with Y)
  a <- as.vector(crossprod(X, Y))
  b <- as.vector(crossprod(M, Y))
  # Normalize initial a and b
  if (sum((X %*% a)^2) == 0) {
    a[1] <- 1  # if zero vector (unlikely), set a trivial non-zero to satisfy norm
  }
  if (sum((M %*% b)^2) == 0) {
    b[1] <- 1
  }
  a <- a / sqrt(sum((X %*% a)^2))
  b <- b / sqrt(sum((M %*% b)^2))
  
  # Initialize auxiliary variables for ADMM
  z_a <- a
  z_b <- b
  u_a <- numeric(m)
  u_b <- numeric(q)
  rho <- 1  # ADMM penalty parameter
  
  # Precompute constants
  P_const <- as.vector(crossprod(X, Y))   # X^T Y (m-vector)
  P2_const <- as.vector(crossprod(M, Y))  # M^T Y (q-vector)
  
  # ADMM iterations
  for (iter in 1:max_iter) {
    a_old <- a
    b_old <- b
    
    ## Update a (optimize a with b fixed)
    inner_iter_a <- 5
    for (j in 1:inner_iter_a) {
      X_a <- as.vector(X %*% a)
      M_b <- as.vector(M %*% b)
      r_yx <- as.numeric(crossprod(a, P_const))
      r_yz <- as.numeric(crossprod(b, P2_const))
      alpha_val <- as.numeric(crossprod(X_a, M_b))
      D_val <- 1 - alpha_val^2
      T_val <- sum(Y^2)
      N_val <- r_yx^2 + r_yz^2 - 2 * alpha_val * r_yx * r_yz
      P_vec <- P_const
      Q_vec <- as.vector(crossprod(X, M_b))
      
      g1 <- -2 * r_yx * P_vec
      w_med <- as.vector(crossprod(M, X_a))
      g2 <- -2 * as.vector(crossprod(X, M %*% w_med))
      dN_da <- 2 * r_yx * P_vec - 2 * r_yz * (r_yx * Q_vec + alpha_val * P_vec)
      if (D_val == 0) D_val <- 1e-6
      g3 <- (2 * T_val * alpha_val * Q_vec) / (D_val^2) - dN_da / (D_val^2) - (4 * alpha_val * N_val * Q_vec) / (D_val^3)
      num_prime <- r_yz * (1 + alpha_val^2) * Q_vec - 2 * alpha_val * r_yx * Q_vec - 2 * alpha_val^2 * (1 - alpha_val^2) * P_vec
      g4 <- - lambda_n * num_prime / (D_val^2)
      g5 <- rho * (a - z_a + u_a)
      
      grad_a <- g1 + g2 + g3 + g4 + g5
      step_size <- 0.1
      a <- a - step_size * grad_a
      if (sum((X %*% a)^2) != 0) {
        a <- a / sqrt(sum((X %*% a)^2))
      }
    }
    
    ## Update b (optimize b with a fixed)
    inner_iter_b <- 5
    for (j in 1:inner_iter_b) {
      X_a <- as.vector(X %*% a)
      M_b <- as.vector(M %*% b)
      r_yx <- as.numeric(crossprod(a, P_const))
      r_yz <- as.numeric(crossprod(b, P2_const))
      alpha_val <- as.numeric(crossprod(X_a, M_b))
      D_val <- 1 - alpha_val^2
      T_val <- sum(Y^2)
      N_val <- r_yx^2 + r_yz^2 - 2 * alpha_val * r_yx * r_yz
      P2_vec <- P2_const
      Q2_vec <- as.vector(crossprod(M, X_a))
      
      dN_db <- 2 * r_yz * P2_vec - 2 * r_yx * (r_yz * Q2_vec + alpha_val * P2_vec)
      if (D_val == 0) D_val <- 1e-6
      g3_b <- (2 * T_val * alpha_val * Q2_vec) / (D_val^2) - dN_db / (D_val^2) - (4 * alpha_val * N_val * Q2_vec) / (D_val^3)
      num_prime_b <- D_val * alpha_val * P2_vec + r_yz * (D_val + 2 * alpha_val^2) * Q2_vec - 2 * r_yx * alpha_val * (D_val + alpha_val^2) * Q2_vec
      g4_b <- - lambda_n * num_prime_b / (D_val^2)
      g5_b <- rho * (b - z_b + u_b)
      
      grad_b <- g3_b + g4_b + g5_b
      step_size <- 0.1
      b <- b - step_size * grad_b
      if (sum((M %*% b)^2) != 0) {
        b <- b / sqrt(sum((M %*% b)^2))
      }
    }
    
    ## Update auxiliary variables z_a, z_b (soft-thresholding for L1)
    soft_threshold <- function(x, thr) { sign(x) * pmax(abs(x) - thr, 0) }
    z_a <- soft_threshold(a + u_a, lambda_a / rho)
    z_b <- soft_threshold(b + u_b, lambda_b / rho)
    
    ## Update dual variables u_a, u_b
    u_a <- u_a + a - z_a
    u_b <- u_b + b - z_b
    
    # Check convergence
    if (sqrt(sum((a - a_old)^2) + sum((b - b_old)^2)) < tol) {
      break
    }
  }
  
  if (sum((X %*% a)^2) != 0) a <- a / sqrt(sum((X %*% a)^2))
  if (sum((M %*% b)^2) != 0) b <- b / sqrt(sum((M %*% b)^2))
  
  return(list(a = a, b = b))
}

###############################################################################
# 3) Cross-Validation for lambda_a, lambda_b
###############################################################################
cv_select_lambda <- function(X, M, Y, lambda_n, lambda_a_seq, lambda_b_seq, K=5) {
  n <- nrow(X)
  folds <- split(sample(1:n), rep(1:K, length.out=n))
  mean_SSR <- matrix(NA, nrow=length(lambda_a_seq), ncol=length(lambda_b_seq))
  dimnames(mean_SSR) <- list(paste0("la=", lambda_a_seq), paste0("lb=", lambda_b_seq))
  
  for (ia in seq_along(lambda_a_seq)) {
    for (ib in seq_along(lambda_b_seq)) {
      lam_a <- lambda_a_seq[ia]
      lam_b <- lambda_b_seq[ib]
      cv_SSR <- numeric(K)
      for (k in 1:K) {
        test_idx <- folds[[k]]
        train_idx <- setdiff(1:n, test_idx)
        X_train <- X[train_idx, , drop=FALSE]
        M_train <- M[train_idx, , drop=FALSE]
        Y_train <- Y[train_idx]
        X_test <- X[test_idx, , drop=FALSE]
        M_test <- M[test_idx, , drop=FALSE]
        Y_test <- Y[test_idx]
        
        model <- optimize_weights(X_train, M_train, Y_train,
                                  lambda_n=lambda_n, 
                                  lambda_a=lam_a, lambda_b=lam_b,
                                  max_iter=50, tol=1e-3)
        a_hat <- model$a
        b_hat <- model$b
        
        x_test <- as.vector(X_test %*% a_hat)
        z_test <- as.vector(M_test %*% b_hat)
        test_fit <- lm(Y_test ~ x_test + z_test)
        rss <- sum(test_fit$residuals^2)
        cv_SSR[k] <- rss
      }
      mean_SSR[ia, ib] <- mean(cv_SSR)
    }
  }
  min_idx <- which(mean_SSR == min(mean_SSR, na.rm=TRUE), arr.ind=TRUE)[1, ]
  best_lambda_a <- lambda_a_seq[min_idx[1]]
  best_lambda_b <- lambda_b_seq[min_idx[2]]
  return(list(best_lambda_a = best_lambda_a,
              best_lambda_b = best_lambda_b,
              cv_error_grid = mean_SSR))
}

###############################################################################
# 4) Simulation wrapper that uses cross-validation and calculates performance
###############################################################################
run_simulation_with_cv <- function(n_runs=10, 
                                   n1=100, m=20, q=10, q_a=3, r=3,
                                   rhoX=0.3, rhoM=0.3,
                                   sigma1=1, sigma2=1,
                                   pathway=c("partial","complete","none"),
                                   lambda_n=1,
                                   lambda_a_seq=c(0, 0.1, 0.2, 0.5),
                                   lambda_b_seq=c(0, 0.1, 0.2, 0.5),
                                   K=5, threshold=1e-3) {
  pathway <- match.arg(pathway)
  results <- data.frame(MP=numeric(0), Accuracy=numeric(0), Precision=numeric(0),
                        Recall=numeric(0), F1=numeric(0))
  
  for (run in 1:n_runs) {
    data <- simulate_data(n1, m, q, q_a, r,
                          rhoX=rhoX, rhoM=rhoM,
                          sigma1=sigma1, sigma2=sigma2,
                          pathway=pathway)
    X <- data$X
    M <- data$M
    Y <- data$Y
    true_active_med <- data$active_mediators
    
    cv_result <- cv_select_lambda(X, M, Y, lambda_n,
                                  lambda_a_seq, lambda_b_seq, K=K)
    best_la <- cv_result$best_lambda_a
    best_lb <- cv_result$best_lambda_b
    
    model <- optimize_weights(X, M, Y,
                              lambda_n=lambda_n,
                              lambda_a=best_la, lambda_b=best_lb,
                              max_iter=50, tol=1e-3)
    a_hat <- model$a
    b_hat <- model$b
    
    # Summarize final model
    x_agg <- as.vector(X %*% a_hat)
    z_agg <- as.vector(M %*% b_hat)
    fit <- lm(Y ~ x_agg + z_agg)
    coef_vals <- coef(fit)
    gamma_hat <- as.numeric(coef_vals["x_agg"])
    eta_hat <- as.numeric(coef_vals["z_agg"])
    alpha_hat <- as.numeric(crossprod(x_agg, z_agg))
    denom <- alpha_hat * eta_hat + gamma_hat
    MP_hat <- if (abs(denom) < 1e-8) 0 else (alpha_hat * eta_hat / denom)
    # Ensure MP is between 0 and 1
    MP_hat <- max(0, min(1, MP_hat))
    
    # Identify active mediators
    pred_active_med <- which(abs(b_hat) > threshold)
    true_active_set <- if (length(true_active_med) > 0) true_active_med else integer(0)
    
    TP <- length(intersect(pred_active_med, true_active_set))
    FP <- length(setdiff(pred_active_med, true_active_set))
    FN <- length(setdiff(true_active_set, pred_active_med))
    TN <- q - (TP + FP + FN)
    
    if ((TP + FP + FN) == 0) {
      # Edge case: no mediators flagged
      accuracy <- 1
      precision <- 1
      recall <- 1
      F1 <- 1
    } else {
      accuracy <- (TP + TN) / q
      precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
      recall <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      if (is.na(precision) || is.na(recall)) {
        F1 <- 0
      } else if ((precision + recall) == 0) {
        F1 <- 0
      } else {
        F1 <- 2 * precision * recall / (precision + recall)
      }
    }
    
    results <- rbind(results,
                     data.frame(MP=MP_hat,
                                Accuracy=accuracy,
                                Precision=precision,
                                Recall=recall,
                                F1=F1))
  }
  
  avg_metrics <- colMeans(results[, -which(colnames(results) == "MP")], na.rm=TRUE)
  avg_MP <- mean(results$MP)
  
  cat(sprintf("Avg MP: %.3f\n", avg_MP))
  cat(sprintf("Avg Accuracy: %.3f\n", avg_metrics["Accuracy"]))
  cat(sprintf("Avg Precision: %.3f\n", avg_metrics["Precision"]))
  cat(sprintf("Avg Recall: %.3f\n", avg_metrics["Recall"]))
  cat(sprintf("Avg F1-score: %.3f\n", avg_metrics["F1"]))
  
  return(results)
}

###############################################################################
# 5) Run multiple scenarios
###############################################################################

set.seed(123)

# Expand your scenario grid to include sigma1 and sigma2:
param_combos <- expand.grid(
  rhoX = c(0, 0.3),
  rhoM = c(0, 0.3),
  sigma1 = c(1, 3),
  sigma2 = c(1, 3),
  mediation_type = c("complete", "partial"),
  stringsAsFactors = FALSE
)

scenario_results <- data.frame()

for (i in seq_len(nrow(param_combos))) {
  px <- param_combos$rhoX[i]
  pm <- param_combos$rhoM[i]
  s1 <- param_combos$sigma1[i]
  s2 <- param_combos$sigma2[i]
  med_type <- param_combos$mediation_type[i]
  
  cat(sprintf("\n=== Running scenario %d of %d ===\n", i, nrow(param_combos)))
  cat(sprintf("   rhoX = %.2f, rhoM = %.2f, sigma1 = %.1f, sigma2 = %.1f, mediation_type = %s\n",
              px, pm, s1, s2, med_type))
  
  # Fix sample size and dimension or set them as you like:
  n1_val <- 200
  m_val <- 20
  q_val <- 20
  
  # Optionally adjust q_a and r based on m and q
  q_a_val <- max(1, round(q_val * 0.25))  # 25% of q
  r_val   <- max(1, round(m_val * 0.25))  # 25% of m
  
  set.seed(123 + i)
  results_df <- run_simulation_with_cv(
    n_runs = 100,
    n1 = n1_val,
    m = m_val,
    q = q_val,
    q_a = q_a_val,
    r = r_val,
    rhoX = px,
    rhoM = pm,
    sigma1 = s1,   # pass sigma1
    sigma2 = s2,   # pass sigma2
    pathway = med_type,
    lambda_n = 5,
    lambda_a_seq = c(0, 0.1, 0.25, 0.5, 1),
    lambda_b_seq = c(0, 0.1, 0.25, 0.5, 1),
    K = 5,
    threshold = 0.006
  )
  
  # Compute summary statistics and, if needed, bias/SD for MP.
  scenario_means <- colMeans(results_df)
  
  # Determine the true mediation proportion for the scenario
  if (med_type == "complete") {
    true_MP <- 1
  } else if (med_type == "partial") {
    true_MP <- 0.5
  }
  
  abs_bias <- abs(results_df$MP - true_MP)
  avg_abs_bias <- mean(abs_bias)
  sd_MP <- sd(results_df$MP)
  
  scenario_results <- rbind(
    scenario_results,
    data.frame(
      rhoX = px,
      rhoM = pm,
      sigma1 = s1,
      sigma2 = s2,
      pathway = med_type,
      MP = scenario_means["MP"],
      Accuracy = scenario_means["Accuracy"],
      Precision = scenario_means["Precision"],
      Recall = scenario_means["Recall"],
      F1 = scenario_means["F1"],
      Abs_Bias = avg_abs_bias,
      SD_MP = sd_MP,
      stringsAsFactors = FALSE
    )
  )
}

cat("\nFinal scenario results:\n")
print(scenario_results)  ## for the noise term 




run_simulation_with_cv <- function(n_runs=10,  
                                   n1, m, q, q_a, r,
                                   rhoX, rhoM, sigma1, sigma2, pathway,
                                   lambda_n, lambda_a_seq, lambda_b_seq,
                                   K, threshold) {
  results <- data.frame(
    MP        = numeric(n_runs),
    Accuracy  = numeric(n_runs),
    Precision = numeric(n_runs),
    Recall    = numeric(n_runs),
    F1        = numeric(n_runs),
    MSE       = numeric(n_runs)   # <--- new column
  )
  
  for(run in seq_len(n_runs)) {
    ## (a) simulate and fit exactly as before…
    data <- simulate_data(n1,m,q,q_a,r,
                          rhoX,rhoM,sigma1,sigma2,pathway)
    X <- data$X; M <- data$M; Y <- data$Y
    cvr <- cv_select_lambda(X,M,Y,lambda_n, lambda_a_seq, lambda_b_seq, K)
    fit_ab <- optimize_weights(X,M,Y,lambda_n,
                               cvr$best_lambda_a, cvr$best_lambda_b)
    a_hat <- fit_ab$a; b_hat <- fit_ab$b
    
    ## (b) fit the aggregate‐mediator model and compute MP
    x_agg <- as.vector(X %*% a_hat)
    z_agg <- as.vector(M %*% b_hat)
    fit    <- lm(Y ~ x_agg + z_agg)
    
    ## ---- here: compute Yhat and its MSE ----
    Y_pred <- predict(fit)                      
    mse_val <- mean((Y - Y_pred)^2)              
    
    ## (c) compute MP, Accuracy, Precision, Recall, F1 exactly as before…
    coefv <- coef(fit)
    gamma_hat <- coefv["x_agg"]
    eta_hat   <- coefv["z_agg"]
    alpha_hat <- as.numeric(crossprod(x_agg, z_agg))
    denom     <- alpha_hat * eta_hat + gamma_hat
    MP_hat    <- if(abs(denom)<1e-8) 0 else (alpha_hat*eta_hat/denom)
    MP_hat    <- max(0, min(1, MP_hat))
    
    pred_active <- which(abs(b_hat) > threshold)
    true_active <- data$active_mediators
    
    TP <- length(intersect(pred_active, true_active))
    FP <- length(setdiff(pred_active, true_active))
    FN <- length(setdiff(true_active, pred_active))
    TN <- q - (TP+FP+FN)
    
    if((TP+FP+FN)==0){
      acc <- prec <- rec <- F1 <- 1
    } else {
      acc  <- (TP+TN)/q
      prec <- if((TP+FP)>0) TP/(TP+FP) else 0
      rec  <- if((TP+FN)>0) TP/(TP+FN) else 0
      F1   <- if((prec+rec)>0) 2*prec*rec/(prec+rec) else 0
    }
    
    ## (d) store everything, **including** mse_val
    results[run, ] <- c(MP_hat, acc, prec, rec, F1, mse_val)
  }
  
  ## Optional: print average MSE too
  cat(sprintf("Average MSE over %d runs: %.4f\n", n_runs, mean(results$MSE)))
  
  return(results)
}

###############################################################################
# (5) Summarize across replicates – also pull in mean_MSE & sd_MSE
###############################################################################
# Expand your scenario grid to include sigma1 and sigma2:
param_combos <- expand.grid(
  rhoX            = c(0, 0.3),
  rhoM            = c(0, 0.3),
  sigma1          = c(1, 3),
  sigma2          = c(1, 3),
  mediation_type  = c("partial"),
  stringsAsFactors = FALSE
)

scenario_results <- data.frame()

for (i in seq_len(nrow(param_combos))) {
  # ←—— HERE IS THE FIX
  px       <- param_combos$rhoX[i]
  pm       <- param_combos$rhoM[i]
  s1       <- param_combos$sigma1[i]
  s2       <- param_combos$sigma2[i]
  med_type <- param_combos$mediation_type[i]
  
  cat(sprintf(
    "\n=== Scenario %d/%d — rhoX=%.1f, rhoM=%.1f, σ1=%.0f, σ2=%.0f, %s\n",
    i, nrow(param_combos), px, pm, s1, s2, med_type
  ))
  
  results_df <- run_simulation_with_cv(
    n_runs       = 20,
    n1           = 200,
    m            = 20,
    q            = 20,
    q_a          = 5,
    r            = 5,
    rhoX         = px,
    rhoM         = pm,
    sigma1       = s1,
    sigma2       = s2,
    pathway      = med_type,
    lambda_n     = 5,
    lambda_a_seq = c(0, 0.1, 0.25, 0.5, 1),
    lambda_b_seq = c(0, 0.1, 0.25, 0.5, 1),
    K            = 5,
    threshold    = 0.006
  )
  
  scenario_results <- rbind(
    scenario_results,
    data.frame(
      rhoX      = px,
      rhoM      = pm,
      sigma1    = s1,
      sigma2    = s2,
      pathway   = med_type,
      MP        = mean(results_df$MP),
      Accuracy  = mean(results_df$Accuracy),
      Precision = mean(results_df$Precision),
      Recall    = mean(results_df$Recall),
      F1        = mean(results_df$F1),
      Mean_MSE  = mean(results_df$MSE),
      SD_MSE    = sd(results_df$MSE),
      stringsAsFactors = FALSE
    )
  )
}

print(scenario_results)


### misspecified model for a 
###############################################################################
#  Miss-specified data–generation function
#  – latent directions 𝐚 and 𝐛 no longer align with the “true” signal       –
#    • X → M loadings α_{ij} are small random values on *every* entry         –
#    • direct effects γ_j and mediator effects η_k are likewise random small  –
###############################################################################
# ---------------------------------------------------------------
#  Simulate high-dimensional mediation data
#  – choose whether α-coefficients (“a’s”) or η-coefficients (“b’s”)
#    are forced to zero ⇢  evaluate robustness
# ---------------------------------------------------------------
## ================================================================
###############################################################################
## 1 ▸  Simulator that toggles α– or η–irrelevance  (unchanged from before) ###
###############################################################################
simulate_data_irrelevant <- function(n1, m, q, r, q_a,
                                     rhoX   = 0,
                                     rhoM   = 0,
                                     sigma1 = 1,
                                     sigma2 = 1,
                                     irrelev = c("none", "a_irrelevant", "b_irrelevant")) {
  
  irrelev <- match.arg(irrelev)
  
  SigmaX <- (1 - rhoX) * diag(m) + rhoX * matrix(1, m, m)
  X      <- unclass(scale(MASS::mvrnorm(n1, rep(0, m), SigmaX)))
  
  SigmaM <- (1 - rhoM) * diag(q) + rhoM * matrix(1, q, q)
  E_M    <- MASS::mvrnorm(n1, rep(0, q), sigma1^2 * SigmaM)
  
  alpha_mat <- matrix(0, m, q)
  eta_vec   <- rep(0, q)
  gamma_vec <- rep(0, m)                # ← zero-vector, **not** scalar 0
  
  for (i in seq_len(r)) if (i <= q_a) {
    alpha_mat[i, i] <- 1
    eta_vec[i]      <- 1
  }
  if (irrelev == "a_irrelevant") alpha_mat[,] <- 0.01
  if (irrelev == "b_irrelevant") eta_vec[]    <- 0.01
  
  M <- unclass(scale(X %*% alpha_mat + E_M))
  
  ## -------- corrected line -----------------------------------------------
  Y <- scale(drop(X %*% gamma_vec + M %*% eta_vec + rnorm(n1, sd = sigma2)))[, 1]
  ## -----------------------------------------------------------------------
  
  list(X = X, M = M, Y = Y,
       alpha_true = alpha_mat, gamma_true = gamma_vec, eta_true = eta_vec)
}
################################################################################
## 2 ▸  Your existing modelling helpers – only the headers shown here.       ##
##     KEEP their internal code exactly as you have it.                       ##
################################################################################



################################################################################
## 3 ▸  run_simulation_with_cv()  – minimal edits marked  ### NEW / CHANGED ###
################################################################################
run_simulation_with_cv <- function(n_runs = 10,
                                   n1, m, q, q_a, r,
                                   rhoX, rhoM,
                                   sigma1, sigma2,
                                   pathway,                  # kept for completeness
                                   irrelev = "none",         # NEW: block-irrelevance
                                   lambda_n,
                                   lambda_a_seq,
                                   lambda_b_seq,
                                   K        = 5,
                                   threshold = 1e-4) {
  
  results <- data.frame(
    MP        = numeric(n_runs),
    Accuracy  = numeric(n_runs),
    Precision = numeric(n_runs),
    Recall    = numeric(n_runs),
    F1        = numeric(n_runs),
    MSE       = numeric(n_runs)
  )
  
  for (run in seq_len(n_runs)) {
    
    ## (a)  simulate one data set   -----------------------------------------
    dat <- simulate_data_irrelevant(n1, m, q, r, q_a,
                                    rhoX, rhoM,
                                    sigma1, sigma2,
                                    irrelev = irrelev)
    
    X <- dat$X         # n × m
    M <- dat$M         # n × q
    Y <- dat$Y         # length n
    true_active <- which(abs(dat$eta_true) > 0)   # ground truth for b
    
    ## (b)  cross-validation to choose λ_a , λ_b  ---------------------------
    best <- cv_select_lambda(X, M, Y,
                             lambda_n,
                             lambda_a_seq,
                             lambda_b_seq,
                             K = K)
    
    lambda_a_best <- best$best_lambda_a   ## <-- adjust if names differ
    lambda_b_best <- best$best_lambda_b   ## <-- adjust if names differ
    
    ## (c)  final fit with chosen lambdas  ----------------------------------
    fit_ab <- optimize_weights(X, M, Y,
                               lambda_n = lambda_n,
                               lambda_a = lambda_a_best,
                               lambda_b = lambda_b_best)
    
    a_hat <- fit_ab$a
    b_hat <- fit_ab$b
    
    ## (d)  aggregate mediators and refit simple Y ~ X* + M*  ---------------
    x_agg <- as.vector(X %*% a_hat)
    z_agg <- as.vector(M %*% b_hat)
    
    lm_fit <- lm(Y ~ x_agg + z_agg)
    Y_pred <- predict(lm_fit)
    mse_val <- mean((Y - Y_pred)^2)
    
    ## mediation-proportion estimate
    coefs <- coef(lm_fit)
    gamma_hat  <- coefs["x_agg"]
    eta_hat    <- coefs["z_agg"]
    alpha_hat  <- as.numeric(crossprod(x_agg, z_agg))  # proxy for α·X loading
    denom      <- alpha_hat * eta_hat + gamma_hat
    MP_hat     <- if (abs(denom) < 1e-8) 0 else (alpha_hat * eta_hat / denom)
    MP_hat     <- max(0, min(1, MP_hat))
    
    ## (e)  variable-selection metrics for η-loadings -----------------------
    pred_active <- which(abs(b_hat) > threshold)
    
    TP <- length(intersect(pred_active, true_active))
    FP <- length(setdiff(pred_active,  true_active))
    FN <- length(setdiff(true_active,  pred_active))
    TN <- q - (TP + FP + FN)
    
    acc  <- if ((TP + FP + FN) == 0) 1 else (TP + TN) / q
    prec <- if ((TP + FP)     > 0) TP / (TP + FP) else 0
    rec  <- if ((TP + FN)     > 0) TP / (TP + FN) else 0
    F1   <- if ((prec + rec)  > 0) 2 * prec * rec / (prec + rec) else 0
    
    ## (f)  stash results ---------------------------------------------------
    results[run, ] <- c(MP_hat, acc, prec, rec, F1, mse_val)
  }
  
  results
}

###############################################################################
## 4 ▸  Full design grid + master loop – identical to previous answer,     ####
##     except we now pass 'irrelev' through.                                 ###
###############################################################################

param_combos <- expand.grid(
  rhoX           = c(0, 0.3),
  rhoM           = c(0, 0.3),
  sigma1         = c(1, 3),
  sigma2         = c(1, 3),
  mediation_type = c("complete"),
  irrelev        = c("a_irrelevant", "b_irrelevant"),
  stringsAsFactors = FALSE
)

scenario_results1 <- data.frame()

for (i in seq_len(nrow(param_combos))) {
  params <- param_combos[i, ]        # pull the row once
  
  cat(sprintf("\n=== Scenario %d/%d — ρX=%.1f, ρM=%.1f, σ1=%.0f, σ2=%.0f, %s, %s\n",
              i, nrow(param_combos),
              params$rhoX, params$rhoM,
              params$sigma1, params$sigma2,
              params$mediation_type, params$irrelev))
  
  res <- run_simulation_with_cv(
    n_runs   = 20,
    n1       = 200, m = 20, q = 20, q_a = 5, r = 5,
    rhoX     = params$rhoX,
    rhoM     = params$rhoM,
    sigma1   = params$sigma1,
    sigma2   = params$sigma2,
    pathway  = params$mediation_type,
    irrelev  = params$irrelev,
    lambda_n = 5,
    lambda_a_seq = c(0, 0.1, 0.25, 0.5, 1),
    lambda_b_seq = c(0, 0.1, 0.25, 0.5, 1),
    K = 5,
    threshold = 0.006)
  
  # build one tidy row
  row_out <- data.frame(t(colMeans(res)),
                        SD_MSE = sd(res$MSE),
                        rhoX   = params$rhoX,
                        rhoM   = params$rhoM,
                        sigma1 = params$sigma1,
                        sigma2 = params$sigma2,
                        pathway = params$mediation_type,
                        irrelev = params$irrelev)
  
  scenario_results1 <- rbind(scenario_results1, row_out)
}

print(scenario_results1)
