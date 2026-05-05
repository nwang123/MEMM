###############################################################################

# Reproduce main-text Table 1: Complete mediation, sigma_Y^2 = 1

###############################################################################

library(MASS)
library(Matrix)
# Optional packages used by comparison methods
# install.packages(c("glmnet", "ncvreg", "SIS", "PMA", "CVXR"))
source("../R/MEMM.R")
set.seed(20260429)

table1_results <- run_table12_summaries(
  nrep = 1000,
  pathways = c("complete"),
  size_grid = data.frame(
    m = c(20L, 50L),
    q = c(20L, 50L),
    stringsAsFactors = FALSE
  ),
  rhoX_values = c(0, 0.3),
  rhoM_values = c(0, 0.3),
  n1 = 200,
  sigma1 = 1,
  sigma2 = 1,
  x_dist = "normal",
  m_noise_dist = "normal",
  df_x = 6,
  df_m = 6,
  c_signal = NULL,
  methods = c(
    "MEMM",
    "SIS_MCP",
    "Pathway",
    "TS",
    "DM1",
    "DM2",
    "SPCMA",
    "FBAS"
  ),
  threshold = 1e-2,
  lambda_fbas = 0.01,
  rho_fbas = 1,
  lambda_n = 1,
  lambda_a_seq = c(0, 0.1, 0.2, 0.5),
  lambda_b_seq = c(0, 0.1, 0.2, 0.5),
  K = 5,
  final_n_restarts = 1,
  proxy_aggregation = "l2",
  truth_mode = "legacy_mediator_only",
  seed = 20260429,
  verbose = TRUE
)

table1_complete <- table1_results$table1_complete
print(table1_complete)

