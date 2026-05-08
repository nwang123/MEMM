###############################################################################

# Reproduce main-text Table 1: Complete mediation, sigma_Y^2 = 1

###############################################################################

library(MASS)
library(Matrix)
# Optional packages used by comparison methods
# install.packages(c("glmnet", "ncvreg", "SIS", "PMA", "CVXR"))
source("../R/MEMM.R")

### Example code for reproducing main-text Table 1.
### This block is provided as an example only and is not run automatically.
### To run it, remove the leading "#" symbols.

# table1_results <- run_table12_summaries(
#
#   ### Number of simulation replications per scenario
#   nrep = 1000,

#   ### Complete mediation only; true MP = 1
#   pathways = c("complete"),
#
#   ### Dimensions considered in Table 1
#   size_grid = data.frame(
#     m = c(20L, 50L),     ### number of exposures
#     q = c(20L, 50L),     ### number of mediators
#     stringsAsFactors = FALSE
#   ),
#
#   ### Exposure and mediator correlation settings
#   rhoX_values = c(0, 0.3),
#   rhoM_values = c(0, 0.3),
#
#   ### Sample size and noise levels
#   n1 = 200,
#   sigma1 = 1,            ### mediator noise scale
#   sigma2 = 1,            ### outcome noise scale; sigma_Y^2 = 1
#
#   ### Distributional settings
#   x_dist = "normal",
#   m_noise_dist = "normal",
#   df_x = 6,              ### unused for normal distribution
#   df_m = 6,              ### unused for normal distribution
#   c_signal = NULL,       ### default signal strength used in main text
#
#   ### Methods included in Table 1
#   methods = c(
#     "MEMM",
#     "SIS_MCP",
#     "Pathway",
#     "TS",
#     "DM1",
#     "DM2",
#     "SPCMA",
#     "FBAS"
#   ),
#
#   ### Support-recovery threshold and comparison-method parameters
#   threshold = 1e-2,
#   lambda_fbas = 0.01,
#   rho_fbas = 1,
#
#   ### MEMM tuning setup
#   lambda_n = 1,
#   lambda_a_seq = c(0, 0.1, 0.2, 0.5),
#   lambda_b_seq = c(0, 0.1, 0.2, 0.5),
#   K = 5,                 ### five-fold cross-validation
#
#   ### Optimization and evaluation settings
#   final_n_restarts = 1,  ### no random restart for main Table 1
#   proxy_aggregation = "l2",
#   truth_mode = "legacy_mediator_only",
#
#   ### Reproducibility
#   seed = 20260429,
#   verbose = TRUE
# )
#
# ### Extract complete mediation results
# table1_complete <- table1_results$table1_complete
#
# ### Print results
# print(table1_complete)
