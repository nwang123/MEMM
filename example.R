# Load required packages
library(MASS)   # For multivariate normal generation
library(glmnet) # For LASSO-like operations if needed

#--- 1. Simulate data --------------------------------------------------------
set.seed(123)

sim_data <- simulate_data(
  n1 = 100,    # Sample size
  m = 10,      # Number of exposures
  q = 8,       # Number of mediators
  q_a = 3,     # Number of active mediators
  r = 2,       # Number of active exposures
  rhoX = 0.3,  # Correlation among exposures
  rhoM = 0.3,  # Correlation among mediators
  sigma1 = 1,  # Noise scale for mediator model
  sigma2 = 1,  # Noise scale for outcome model
  pathway = "partial"
)

# The simulated data includes:
# sim_data$X : exposure matrix (n × m)
# sim_data$M : mediator matrix (n × q)
# sim_data$Y : outcome vector (n × 1)
# sim_data$alpha, sim_data$eta : true coefficient matrices

#--- 2. Cross-validate LASSO penalties ---------------------------------------
cv_results <- cv_select_lambda(
  X = sim_data$X,
  M = sim_data$M,
  Y = sim_data$Y,
  lambda_seq = c(0, 0.1, 0.25, 0.5, 1), # tuning grid
  folds = 5
)

cat("Selected λ_a:", cv_results$best_lambda_a, "\n")
cat("Selected λ_b:", cv_results$best_lambda_b, "\n")

#--- 3. Fit model via ADMM optimization -------------------------------------
fit <- optimize_weights(
  X = sim_data$X,
  M = sim_data$M,
  Y = sim_data$Y,
  lambda_a = cv_results$best_lambda_a,
  lambda_b = cv_results$best_lambda_b,
  tol = 1e-4,
  max_iter = 1000
)

# Estimated projection directions:
a_hat <- fit$a
b_hat <- fit$b

#--- 4. Evaluate performance ------------------------------------------------
metrics <- evaluate_performance(
  a_true = sim_data$a_true,
  b_true = sim_data$b_true,
  a_est = a_hat,
  b_est = b_hat
)

print(metrics)

#--- 5. (Optional) Run the automated wrapper ---------------------------------
summary_results <- run_simulation_with_cv(
  n_runs = 5,
  n1 = 100,
  m = 10,
  q = 8,
  q_a = 3,
  r = 2,
  rhoX = 0.3,
  rhoM = 0.3,
  sigma1 = 1,
  sigma2 = 1,
  pathway = "partial"
)

print(summary_results)
