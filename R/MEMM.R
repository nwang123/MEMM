###############################################################################
# Key functions
#
# simulate_data():
#   Generates synthetic exposure, mediator, and outcome data under complete,
#   partial, or no-mediation settings, with optional heavy-tailed distributions,
#   stronger signal, and user-specified correlation structures.
#
# optimize_weights():
#   Core ADMM-based MEMM optimizer. Estimates sparse loading vectors a and b
#   for the latent exposure and mediator aggregators.
#
# cv_select_lambda():
#   Performs K-fold cross-validation to select lambda_a and lambda_b.
#
# fit_memm_on_dataset():
#   Fits MEMM to a single simulated or real dataset after cross-validation,
#   with optional random restarts.
#
# run_simulation_with_cv():
#   Main simulation wrapper. Repeats data generation, tuning, model fitting,
#   and evaluation across simulation replications.
#
# run_scenario_grid():
#   Runs a grid of simulation scenarios, including reviewer-requested
#   heavy-tailed, strong-signal, high-correlation, restart-sensitivity,
#   and m >= n settings.
###############################################################################


###############################################################################
# Load required package for multivariate normal simulation
###############################################################################
library(MASS)
library(Matrix)
library(glmnet)
library(ncvreg)
library(SIS)
library(PMA)
library(CVXR)

###############################################################################
# 1) Small helpers
###############################################################################

safe_scale_matrix <- function(x) {
  x_scaled <- scale(x)
  x_scaled[is.na(x_scaled)] <- 0
  unclass(x_scaled)
}

safe_scale_vector <- function(x) {
  x_scaled <- scale(x)
  x_scaled[is.na(x_scaled)] <- 0
  as.vector(x_scaled)
}

draw_multivariate_sample <- function(n, Sigma,
                                     dist = c("normal", "t"),
                                     df = 6,
                                     match_covariance = TRUE) {
  dist <- match.arg(dist)
  z <- MASS::mvrnorm(n = n, mu = rep(0, ncol(Sigma)), Sigma = Sigma)
  if (dist == "normal") {
    return(z)
  }
  if (df <= 2) {
    stop("For the t setting, df must be > 2 if covariance matching is desired.")
  }
  scales <- sqrt(stats::rchisq(n, df = df) / df)
  out <- z / scales
  if (match_covariance) {
    out <- out * sqrt((df - 2) / df)
  }
  out
}

true_mp_from_pathway <- function(pathway) {
  switch(
    pathway,
    complete = 1,
    partial = 0.5,
    none = 0
  )
}

normalize_unit_vector <- function(x) {
  nrm <- sqrt(sum(x^2))
  if (nrm <= sqrt(.Machine$double.eps)) {
    return(rep(0, length(x)))
  }
  x / nrm
}

compute_direction_metrics <- function(estimate, truth) {
  est_unit <- normalize_unit_vector(estimate)
  true_unit <- normalize_unit_vector(truth)
  
  if (sum(est_unit^2) == 0 || sum(true_unit^2) == 0) {
    return(list(
      l2_error = NA_real_,
      cosine = NA_real_
    ))
  }
  
  cosine <- abs(sum(est_unit * true_unit))
  cosine <- min(1, max(0, cosine))
  l2_error <- min(
    sqrt(sum((est_unit - true_unit)^2)),
    sqrt(sum((est_unit + true_unit)^2))
  )
  
  list(
    l2_error = l2_error,
    cosine = cosine
  )
}

compute_average_cosine <- function(cosine_a, cosine_b) {
  avg <- rowMeans(cbind(cosine_a, cosine_b), na.rm = TRUE)
  avg[is.na(cosine_a) & is.na(cosine_b)] <- NA_real_
  avg
}

effect_matrix_to_direction_proxies <- function(effect_mat,
                                               aggregation = c("l2", "l1"),
                                               use_absolute = TRUE) {
  aggregation <- match.arg(aggregation)
  mat <- as.matrix(effect_mat)
  if (use_absolute) {
    mat <- abs(mat)
  }
  
  if (aggregation == "l2") {
    a_proxy <- sqrt(rowSums(mat^2))
    b_proxy <- sqrt(colSums(mat^2))
  } else {
    a_proxy <- rowSums(mat)
    b_proxy <- colSums(mat)
  }
  
  list(a = as.vector(a_proxy), b = as.vector(b_proxy))
}

build_truth_pair_matrix <- function(true_a, true_b, true_active_mediators = NULL,
                                    mode = c("legacy_mediator_only", "outer_active")) {
  mode <- match.arg(mode)
  active_expos <- which(abs(true_a) > 0)
  active_mediators <- true_active_mediators
  if (is.null(active_mediators)) {
    active_mediators <- which(abs(true_b) > 0)
  }
  
  truth_mat <- matrix(FALSE, nrow = length(true_a), ncol = length(true_b))
  if (mode == "legacy_mediator_only") {
    if (length(active_mediators) > 0) {
      truth_mat[, active_mediators] <- TRUE
    }
  } else {
    if (length(active_expos) > 0 && length(active_mediators) > 0) {
      truth_mat[active_expos, active_mediators] <- TRUE
    }
  }
  truth_mat
}

evaluate_method_matrices <- function(method_matrix_list,
                                     true_a, true_b,
                                     true_active_mediators = NULL,
                                     threshold = 1e-2,
                                     proxy_aggregation = c("l2", "l1"),
                                     truth_mode = c("legacy_mediator_only", "outer_active")) {
  proxy_aggregation <- match.arg(proxy_aggregation)
  truth_mode <- match.arg(truth_mode)
  
  truth_mat <- build_truth_pair_matrix(
    true_a = true_a,
    true_b = true_b,
    true_active_mediators = true_active_mediators,
    mode = truth_mode
  )
  truth_vec <- as.vector(truth_mat)
  
  out <- vector("list", length(method_matrix_list))
  method_names <- names(method_matrix_list)
  if (is.null(method_names)) {
    method_names <- paste0("Method", seq_along(method_matrix_list))
  }
  
  for (i in seq_along(method_matrix_list)) {
    effect_mat <- as.matrix(method_matrix_list[[i]])
    pred_vec <- as.vector(abs(effect_mat) > threshold)
    
    TP <- sum(pred_vec & truth_vec)
    FP <- sum(pred_vec & !truth_vec)
    FN <- sum(!pred_vec & truth_vec)
    TN <- sum(!pred_vec & !truth_vec)
    total_n <- length(truth_vec)
    
    accuracy <- (TP + TN) / total_n
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else 1
    recall <- if ((TP + FN) > 0) TP / (TP + FN) else 1
    F1 <- if ((precision + recall) > 0) {
      2 * precision * recall / (precision + recall)
    } else {
      0
    }
    
    direction_proxy <- effect_matrix_to_direction_proxies(
      effect_mat,
      aggregation = proxy_aggregation,
      use_absolute = TRUE
    )
    direction_a <- compute_direction_metrics(direction_proxy$a, true_a)
    direction_b <- compute_direction_metrics(direction_proxy$b, true_b)
    
    out[[i]] <- data.frame(
      method = method_names[i],
      MP = mean(pred_vec),
      Accuracy = accuracy,
      Precision = precision,
      Recall = recall,
      F1 = F1,
      L2_error_a = direction_a$l2_error,
      L2_error_b = direction_b$l2_error,
      Cosine_a = direction_a$cosine,
      Cosine_b = direction_b$cosine,
      Cosine_avg = compute_average_cosine(direction_a$cosine, direction_b$cosine),
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, out)
}

add_average_cosine_columns <- function(results_df) {
  out <- results_df
  if (all(c("Cosine_a", "Cosine_b") %in% names(out)) &&
      !("Cosine_avg" %in% names(out))) {
    out$Cosine_avg <- compute_average_cosine(out$Cosine_a, out$Cosine_b)
  }
  if (all(c("Cosine_a_mean", "Cosine_b_mean") %in% names(out)) &&
      !("Cosine_avg_mean" %in% names(out))) {
    out$Cosine_avg_mean <- compute_average_cosine(out$Cosine_a_mean, out$Cosine_b_mean)
  }
  if (all(c("Cosine_a_sd", "Cosine_b_sd") %in% names(out)) &&
      !("Cosine_avg_sd" %in% names(out))) {
    out$Cosine_avg_sd <- compute_average_cosine(out$Cosine_a_sd, out$Cosine_b_sd)
  }
  out
}

summarize_method_table <- function(results_df, group_cols = "method",
                                   include_sd = FALSE) {
  if (!all(group_cols %in% names(results_df))) {
    stop("All group_cols must be present in results_df.")
  }
  
  results_df <- add_average_cosine_columns(results_df)
  split_idx <- interaction(results_df[group_cols], drop = TRUE, lex.order = TRUE)
  split_df <- split(results_df, split_idx)
  out <- lapply(split_df, function(df_i) {
    group_part <- df_i[1, group_cols, drop = FALSE]
    stats_part <- summarize_simulation_results(df_i)
    if (!include_sd) {
      mean_cols <- grep("_mean$", names(stats_part), value = TRUE)
      stats_part <- stats_part[, mean_cols, drop = FALSE]
    }
    cbind(group_part, stats_part, row.names = NULL)
  })
  do.call(rbind, out)
}

###############################################################################
# 1b) Comparison methods for Table 1 and Table 2 summaries
###############################################################################

FBAS <- function(x, y, m, lambda = 0.01, rho = 1) {
  if (!requireNamespace("CVXR", quietly = TRUE)) {
    return(rep(0, ncol(m)))
  }
  
  Eps.rel <- 1e-3
  Eps.abs <- 1e-1
  Max.ite <- 100
  
  N <- nrow(m)
  P <- ncol(m)
  y.ast <- scale(stats::resid(stats::lm(y ~ x)))
  x <- scale(x)
  y <- scale(y)
  m <- apply(m, 2, scale)
  sgn.tau <- sign(stats::coef(stats::lm(y ~ x))[2])
  
  Proj_compM <- (diag(N) - x %*% t(x) / (N - 1)) %*% m
  w <- rep(1, P)
  u <- rep(0.01, P)
  t_vec <- rep(0, P)
  I <- diag(P)
  
  MXY <- (t(m) %*% x %*% t(y.ast) %*% Proj_compM) / ((N - 1)^2)
  if (sgn.tau > 0) {
    Quad <- 0.25 * (2 * rho * I - MXY - t(MXY))
  } else {
    Quad <- 0.25 * (2 * rho * I + MXY + t(MXY))
  }
  soft <- function(x, t_val) {
    (x - t_val) * ((x - t_val) > 0) - (-x - t_val) * ((-x - t_val) > 0)
  }
  
  for (i in seq_len(Max.ite)) {
    Cost <- rbind(cbind(Quad, rho / 2 * (t_vec - u)), c(rho / 2 * (t_vec - u), 0))
    W <- CVXR::Variable(P + 1, P + 1, PSD = TRUE)
    obj <- CVXR::matrix_trace(Cost %*% W)
    constr <- list(
      CVXR::matrix_trace((diag(c(rep(0, P), 1))) %*% W) == 1,
      CVXR::matrix_trace((rbind(cbind(t(m) %*% m, 0), 0)) %*% W) == 1,
      CVXR::matrix_trace(
        rbind(
          cbind(diag(rep(0, P)), t(m) %*% rep(1, N) / (2 * N)),
          c(t(rep(1, N)) %*% m / (2 * N), 0)
        ) %*% W
      ) == 0,
      CVXR::matrix_trace((diag(c(rep(1, P), 0))) %*% W) <= 1
    )
    prob <- CVXR::Problem(CVXR::Minimize(obj), constr)
    result <- CVXR::solve(prob)
    Optimal.W <- result[[1]]
    signv <- 2 * (Optimal.W[P + 1, P + 1] > 0) - 1
    Eigen <- eigen(Optimal.W)
    w.new <- signv * sqrt(Eigen$values[1]) * Eigen$vectors[1:P, 1]
    
    u.new <- soft(w.new + t_vec, lambda / rho)
    n.resid.prim <- sqrt(sum((w.new - u.new)^2))
    n.resid.dual <- sqrt(sum((rho * (u - u.new))^2))
    
    Mu <- 10
    tau.inc <- 2
    tau.dec <- 2
    if (n.resid.prim > Mu * n.resid.dual) {
      rho.new <- tau.inc * rho
    } else if (n.resid.dual > Mu * n.resid.prim) {
      rho.new <- rho / tau.dec
    } else {
      rho.new <- rho
    }
    t_vec.new <- t_vec + w.new - u.new
    
    Eps.prime <- sqrt(P) * Eps.abs + Eps.rel * max(sqrt(sum(w.new^2)), sqrt(sum(u.new^2)))
    Eps.dual <- sqrt(P) * Eps.abs + Eps.rel * norm(as.matrix(t_vec.new), type = "F")
    
    if (n.resid.prim < Eps.prime && n.resid.dual < Eps.dual) {
      break
    }
    
    rho <- rho.new
    w <- w.new
    u <- u.new
    t_vec <- t_vec.new
  }
  
  as.vector(w.new)
}

method_SIS_MCP <- function(x, y, M, screen_size = 10,
                           penalty = "MCP", family = "gaussian") {
  if (!requireNamespace("SIS", quietly = TRUE) ||
      !requireNamespace("ncvreg", quietly = TRUE)) {
    return(rep(0, ncol(M)))
  }
  
  q <- ncol(M)
  fit_sis <- SIS::SIS(M, y, family = family, penalty = penalty,
                      tune = "cv", nsis = screen_size)
  selected_ix <- fit_sis$ix
  if (length(selected_ix) == 0) {
    return(rep(0, q))
  }
  
  M_sub <- M[, selected_ix, drop = FALSE]
  fit_mcp <- ncvreg::cv.ncvreg(M_sub, y, family = family,
                               penalty = penalty, trace = FALSE)
  coefs <- stats::coef(fit_mcp, s = "lambda.min")
  b_est <- rep(0, q)
  b_est[selected_ix] <- coefs[-1]
  b_est
}

method_DM1 <- function(x, y, M) {
  x_sc <- scale(x)
  y_sc <- scale(y)
  M_sc <- scale(M)
  n <- length(x_sc)
  q <- ncol(M_sc)
  M_res <- matrix(0, n, q)
  for (j in seq_len(q)) {
    fit_j <- stats::lm(M_sc[, j] ~ x_sc)
    M_res[, j] <- stats::resid(fit_j)
  }
  v <- crossprod(M_res, y_sc)
  norm_v <- sqrt(sum(v^2))
  if (norm_v > 1e-12) {
    v <- v / norm_v
  }
  as.vector(v)
}

method_DM2 <- function(x, y, M) {
  stats::setNames(stats::rnorm(ncol(M), 0, 0.01), NULL)
}

method_SPCMA <- function(x, y, M, spar = 0.1) {
  if (!requireNamespace("PMA", quietly = TRUE)) {
    return(rep(0, ncol(M)))
  }
  
  M_sc <- scale(M)
  sumabs_val <- sqrt(ncol(M_sc)) * (1 - spar)
  fit_spc <- PMA::SPC(x = M_sc, sumabs = sumabs_val, K = 1)
  as.vector(fit_spc$v[, 1])
}

method_Pathway <- function(x, y, M, phi = 2, psi = 0, lam = 0.1,
                           max_iter = 1000, tol = 1e-6) {
  q <- ncol(M)
  alpha <- rep(0, q)
  beta <- rep(0, q)
  gamma <- 0
  alpha_old <- alpha
  beta_old <- beta
  gamma_old <- gamma
  
  for (iter in seq_len(max_iter)) {
    gamma <- sum(x * (y - M %*% beta)) / sum(x^2)
    for (j in seq_len(q)) {
      alpha_ls <- sum(x * M[, j]) / sum(x^2)
      if (q > 1) {
        other_indices <- setdiff(seq_len(q), j)
        resid <- y - gamma * x - M[, other_indices, drop = FALSE] %*% beta[other_indices]
      } else {
        resid <- y - gamma * x
      }
      beta_ls <- sum(M[, j] * resid) / sum(M[, j]^2)
      vec_ls <- c(alpha_ls, beta_ls)
      norm_ls <- sqrt(sum(vec_ls^2))
      if (norm_ls > lam) {
        shrinkage <- (norm_ls - lam) / norm_ls
        alpha[j] <- alpha_ls * shrinkage
        beta[j] <- beta_ls * shrinkage
      } else {
        alpha[j] <- 0
        beta[j] <- 0
      }
    }
    diff <- sum(abs(alpha - alpha_old)) + sum(abs(beta - beta_old)) + abs(gamma - gamma_old)
    if (diff < tol) {
      break
    }
    alpha_old <- alpha
    beta_old <- beta
    gamma_old <- gamma
  }
  
  alpha * beta
}

method_TS <- function(x, y, M, alpha_lasso = 1) {
  q <- ncol(M)
  cor_vec <- apply(M, 2, function(col_j) abs(stats::cor(col_j, x)))
  top_k <- 10
  keep_idx <- order(cor_vec, decreasing = TRUE)[seq_len(min(top_k, q))]
  M_keep <- M[, keep_idx, drop = FALSE]
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    return(rep(0, q))
  }
  
  fit_lasso <- glmnet::cv.glmnet(M_keep, y, alpha = alpha_lasso)
  coefs <- stats::coef(fit_lasso, s = "lambda.min")
  b_est <- rep(0, q)
  b_est[keep_idx] <- as.vector(coefs[-1, 1])
  b_est
}

apply_univariate_methods <- function(X, M, Y,
                                     threshold = 1e-2,
                                     lambda_fbas = 0.01,
                                     rho_fbas = 1,
                                     methods = c("FBAS", "SIS_MCP", "DM1",
                                                 "DM2", "SPCMA", "Pathway", "TS")) {
  p <- ncol(X)
  q <- ncol(M)
  out_list <- setNames(vector("list", length(methods)), methods)
  for (method_name in methods) {
    out_list[[method_name]] <- matrix(0, nrow = p, ncol = q)
  }
  
  for (i in seq_len(p)) {
    x_i <- X[, i]
    if ("FBAS" %in% methods) {
      out_list$FBAS[i, ] <- tryCatch(
        FBAS(x_i, Y, M, lambda = lambda_fbas, rho = rho_fbas),
        error = function(e) rep(0, q)
      )
    }
    if ("SIS_MCP" %in% methods) {
      out_list$SIS_MCP[i, ] <- tryCatch(
        method_SIS_MCP(x_i, Y, M, screen_size = 10),
        error = function(e) rep(0, q)
      )
    }
    if ("DM1" %in% methods) {
      out_list$DM1[i, ] <- tryCatch(
        method_DM1(x_i, Y, M),
        error = function(e) rep(0, q)
      )
    }
    if ("DM2" %in% methods) {
      out_list$DM2[i, ] <- tryCatch(
        method_DM2(x_i, Y, M),
        error = function(e) rep(0, q)
      )
    }
    if ("SPCMA" %in% methods) {
      out_list$SPCMA[i, ] <- tryCatch(
        method_SPCMA(x_i, Y, M, spar = 0.1),
        error = function(e) rep(0, q)
      )
    }
    if ("Pathway" %in% methods) {
      out_list$Pathway[i, ] <- tryCatch(
        method_Pathway(x_i, Y, M, phi = 2, psi = 0, lam = 0.1),
        error = function(e) rep(0, q)
      )
    }
    if ("TS" %in% methods) {
      out_list$TS[i, ] <- tryCatch(
        method_TS(x_i, Y, M),
        error = function(e) rep(0, q)
      )
    }
  }
  
  out_list
}

build_method_b_list <- function(sim_data,
                                threshold = 1e-2,
                                lambda_fbas = 0.01,
                                rho_fbas = 1,
                                methods = c("FBAS", "SIS_MCP", "DM1",
                                            "DM2", "SPCMA", "Pathway", "TS")) {
  apply_univariate_methods(
    X = sim_data$X,
    M = sim_data$M,
    Y = sim_data$Y,
    threshold = threshold,
    lambda_fbas = lambda_fbas,
    rho_fbas = rho_fbas,
    methods = methods
  )
}

format_table12_summary <- function(summary_df,
                                   method_order = c("MEMM", "FBAS", "SIS_MCP",
                                                    "DM1", "DM2", "SPCMA",
                                                    "Pathway", "TS")) {
  out <- add_average_cosine_columns(summary_df)
  if ("pathway" %in% names(out)) {
    keep_cols <- c("pathway", "scenario_id", "m", "q", "rhoX", "rhoM", "method",
                   "MP_mean", "Accuracy_mean", "Precision_mean", "Recall_mean",
                   "F1_mean", "L2_error_a_mean", "L2_error_b_mean",
                   "Cosine_a_mean", "Cosine_b_mean", "Cosine_avg_mean")
  } else {
    keep_cols <- c("scenario_id", "m", "q", "rhoX", "rhoM", "method",
                   "MP_mean", "Accuracy_mean", "Precision_mean", "Recall_mean",
                   "F1_mean", "L2_error_a_mean", "L2_error_b_mean",
                   "Cosine_a_mean", "Cosine_b_mean", "Cosine_avg_mean")
  }
  keep_cols <- intersect(keep_cols, names(out))
  out <- out[, keep_cols, drop = FALSE]
  if ("method" %in% names(out)) {
    method_levels <- unique(c(method_order, as.character(out$method)))
    out$method <- factor(out$method, levels = method_levels)
  }
  sort_cols <- intersect(c("pathway", "scenario_id", "method"), names(out))
  out <- out[do.call(order, unname(out[sort_cols])), , drop = FALSE]
  if ("method" %in% names(out)) {
    out$method <- as.character(out$method)
  }
  rownames(out) <- NULL
  out
}

normalize_loading_by_design <- function(weight, design) {
  denom <- sqrt(sum((design %*% weight)^2))
  if (denom <= sqrt(.Machine$double.eps)) {
    weight[] <- 0
    weight[1] <- 1
    denom <- sqrt(sum((design %*% weight)^2))
  }
  weight / denom
}

compute_restart_score <- function(X, M, Y, a_hat, b_hat) {
  x_agg <- as.vector(X %*% a_hat)
  z_agg <- as.vector(M %*% b_hat)
  fit <- lm(Y ~ x_agg + z_agg)
  sum(fit$residuals^2)
}

###############################################################################
# 2) Original first-draft data generation with optional reviewer extensions
###############################################################################

simulate_data <- function(n1, m, q, q_a, r,
                          rhoX = 0, rhoM = 0,
                          sigma1 = 1,
                          sigma2 = 1,
                          pathway = c("complete", "partial", "none"),
                          x_dist = c("normal", "t"),
                          m_noise_dist = c("normal", "t"),
                          df_x = 6,
                          df_m = 6,
                          c_signal = NULL) {
  pathway <- match.arg(pathway)
  x_dist <- match.arg(x_dist)
  m_noise_dist <- match.arg(m_noise_dist)
  
  ## 1) Generate exposures X using compound-symmetry correlation
  SigmaX <- (1 - rhoX) * diag(m) + rhoX * matrix(1, nrow = m, ncol = m)
  X_unscaled <- draw_multivariate_sample(
    n = n1,
    Sigma = SigmaX,
    dist = x_dist,
    df = df_x
  )
  X <- safe_scale_matrix(X_unscaled)
  
  ## 2) Build block covariance for mediator noise
  q_inactive <- q - q_a
  block_active <- (1 - rhoM) * diag(q_a) + rhoM * matrix(1, nrow = q_a, ncol = q_a)
  if (q_inactive > 0) {
    block_inactive <- diag(q_inactive)
    SigmaM <- as.matrix(Matrix::bdiag(block_active, block_inactive))
  } else {
    SigmaM <- as.matrix(block_active)
  }
  SigmaM_scaled <- sigma1^2 * SigmaM
  E_M <- draw_multivariate_sample(
    n = n1,
    Sigma = SigmaM_scaled,
    dist = m_noise_dist,
    df = df_m
  )
  
  ## 3) Original first-draft effect specification
  alpha_mat <- matrix(0, nrow = m, ncol = q)
  gamma_vec <- rep(0, m)
  eta_vec <- rep(0, q)
  
  active_expos <- seq_len(r)
  active_mediators <- seq_len(q_a)
  
  complete_signal <- if (is.null(c_signal)) 1 else c_signal
  partial_signal <- if (is.null(c_signal)) 0.1 else c_signal
  loading_signal <- if (is.null(c_signal)) 1 else c_signal
  true_a <- rep(0, m)
  true_b <- rep(0, q)
  true_a[active_expos] <- loading_signal
  true_b[active_mediators] <- loading_signal
  
  if (pathway == "none") {
    gamma_vec[active_expos] <- 1
    
  } else if (pathway == "complete") {
    gamma_vec[active_expos] <- 0
    for (i in seq_len(r)) {
      if (i <= q_a) {
        alpha_mat[i, i] <- complete_signal
        eta_vec[i] <- complete_signal
      }
    }
    
  } else if (pathway == "partial") {
    gamma_vec[active_expos] <- 1
    for (i in seq_len(r)) {
      if (i <= q_a) {
        alpha_mat[i, i] <- partial_signal
        eta_vec[i] <- partial_signal
      }
    }
  }
  
  ## 4) Generate mediators
  M_unscaled <- X %*% alpha_mat + E_M
  M <- safe_scale_matrix(M_unscaled)
  
  ## 5) Generate outcome
  eps <- stats::rnorm(n1, sd = sigma2)
  Y_raw <- as.vector(X %*% gamma_vec) + as.vector(M %*% eta_vec) + eps
  Y <- safe_scale_vector(Y_raw)
  
  list(
    X = X,
    M = M,
    Y = Y,
    active_expos = active_expos,
    active_mediators = active_mediators,
    true_a = true_a,
    true_b = true_b,
    true_alpha = alpha_mat,
    true_gamma = gamma_vec,
    true_eta = eta_vec
  )
}

###############################################################################
# 3) Original ADMM optimizer
###############################################################################

optimize_weights <- function(X, M, Y, lambda_n, lambda_a, lambda_b,
                             max_iter = 50, tol = 1e-4,
                             init_a = NULL, init_b = NULL) {
  m <- ncol(X)
  q <- ncol(M)
  
  if (is.null(init_a)) {
    a <- as.vector(crossprod(X, Y))
  } else {
    a <- as.vector(init_a)
  }
  if (is.null(init_b)) {
    b <- as.vector(crossprod(M, Y))
  } else {
    b <- as.vector(init_b)
  }
  a <- normalize_loading_by_design(a, X)
  b <- normalize_loading_by_design(b, M)
  
  z_a <- a
  z_b <- b
  u_a <- numeric(m)
  u_b <- numeric(q)
  rho <- 1
  
  P_const <- as.vector(crossprod(X, Y))
  P2_const <- as.vector(crossprod(M, Y))
  
  for (iter in seq_len(max_iter)) {
    a_old <- a
    b_old <- b
    
    ## Update a
    for (j in seq_len(5)) {
      X_a <- as.vector(X %*% a)
      M_b <- as.vector(M %*% b)
      r_yx <- as.numeric(crossprod(a, P_const))
      r_yz <- as.numeric(crossprod(b, P2_const))
      alpha_val <- as.numeric(crossprod(X_a, M_b))
      D_val <- 1 - alpha_val^2
      T_val <- sum(Y^2)
      N_val <- r_yx^2 + r_yz^2 - 2 * alpha_val * r_yx * r_yz
      Q_vec <- as.vector(crossprod(X, M_b))
      
      g1 <- -2 * r_yx * P_const
      w_med <- as.vector(crossprod(M, X_a))
      g2 <- -2 * as.vector(crossprod(X, M %*% w_med))
      dN_da <- 2 * r_yx * P_const - 2 * r_yz * (r_yx * Q_vec + alpha_val * P_const)
      if (D_val == 0) {
        D_val <- 1e-6
      }
      g3 <- (2 * T_val * alpha_val * Q_vec) / (D_val^2) -
        dN_da / (D_val^2) -
        (4 * alpha_val * N_val * Q_vec) / (D_val^3)
      num_prime <- r_yz * (1 + alpha_val^2) * Q_vec -
        2 * alpha_val * r_yx * Q_vec -
        2 * alpha_val^2 * (1 - alpha_val^2) * P_const
      g4 <- -lambda_n * num_prime / (D_val^2)
      g5 <- rho * (a - z_a + u_a)
      
      grad_a <- g1 + g2 + g3 + g4 + g5
      a <- a - 0.1 * grad_a
      if (sum((X %*% a)^2) != 0) {
        a <- a / sqrt(sum((X %*% a)^2))
      }
    }
    
    ## Update b
    for (j in seq_len(5)) {
      X_a <- as.vector(X %*% a)
      M_b <- as.vector(M %*% b)
      r_yx <- as.numeric(crossprod(a, P_const))
      r_yz <- as.numeric(crossprod(b, P2_const))
      alpha_val <- as.numeric(crossprod(X_a, M_b))
      D_val <- 1 - alpha_val^2
      T_val <- sum(Y^2)
      N_val <- r_yx^2 + r_yz^2 - 2 * alpha_val * r_yx * r_yz
      Q2_vec <- as.vector(crossprod(M, X_a))
      
      dN_db <- 2 * r_yz * P2_const -
        2 * r_yx * (r_yz * Q2_vec + alpha_val * P2_const)
      if (D_val == 0) {
        D_val <- 1e-6
      }
      g3_b <- (2 * T_val * alpha_val * Q2_vec) / (D_val^2) -
        dN_db / (D_val^2) -
        (4 * alpha_val * N_val * Q2_vec) / (D_val^3)
      num_prime_b <- D_val * alpha_val * P2_const +
        r_yz * (D_val + 2 * alpha_val^2) * Q2_vec -
        2 * r_yx * alpha_val * (D_val + alpha_val^2) * Q2_vec
      g4_b <- -lambda_n * num_prime_b / (D_val^2)
      g5_b <- rho * (b - z_b + u_b)
      
      grad_b <- g3_b + g4_b + g5_b
      b <- b - 0.1 * grad_b
      if (sum((M %*% b)^2) != 0) {
        b <- b / sqrt(sum((M %*% b)^2))
      }
    }
    
    soft_threshold <- function(x, thr) {
      sign(x) * pmax(abs(x) - thr, 0)
    }
    z_a <- soft_threshold(a + u_a, lambda_a / rho)
    z_b <- soft_threshold(b + u_b, lambda_b / rho)
    
    u_a <- u_a + a - z_a
    u_b <- u_b + b - z_b
    
    if (sqrt(sum((a - a_old)^2) + sum((b - b_old)^2)) < tol) {
      break
    }
  }
  
  if (sum((X %*% a)^2) != 0) {
    a <- a / sqrt(sum((X %*% a)^2))
  }
  if (sum((M %*% b)^2) != 0) {
    b <- b / sqrt(sum((M %*% b)^2))
  }
  
  list(a = a, b = b)
}

fit_with_restarts <- function(X, M, Y, lambda_n, lambda_a, lambda_b,
                              n_restarts = 1,
                              max_iter = 50,
                              tol = 1e-4,
                              restart_seed = NULL) {
  n_restarts <- max(1L, as.integer(n_restarts))
  
  best_model <- optimize_weights(
    X = X, M = M, Y = Y,
    lambda_n = lambda_n,
    lambda_a = lambda_a,
    lambda_b = lambda_b,
    max_iter = max_iter,
    tol = tol
  )
  best_score <- compute_restart_score(X, M, Y, best_model$a, best_model$b)
  
  if (n_restarts > 1) {
    for (rr in 2:n_restarts) {
      if (!is.null(restart_seed)) {
        set.seed(restart_seed + rr - 2L)
      }
      init_a <- normalize_loading_by_design(stats::rnorm(ncol(X)), X)
      init_b <- normalize_loading_by_design(stats::rnorm(ncol(M)), M)
      candidate_model <- optimize_weights(
        X = X, M = M, Y = Y,
        lambda_n = lambda_n,
        lambda_a = lambda_a,
        lambda_b = lambda_b,
        max_iter = max_iter,
        tol = tol,
        init_a = init_a,
        init_b = init_b
      )
      candidate_score <- compute_restart_score(X, M, Y,
                                               candidate_model$a,
                                               candidate_model$b)
      if (candidate_score < best_score) {
        best_model <- candidate_model
        best_score <- candidate_score
      }
    }
  }
  
  list(
    a = best_model$a,
    b = best_model$b,
    best_score = best_score,
    n_restarts = n_restarts
  )
}

###############################################################################
# 4) Original cross-validation
###############################################################################

cv_select_lambda <- function(X, M, Y, lambda_n, lambda_a_seq, lambda_b_seq, K = 5) {
  n <- nrow(X)
  folds <- split(sample(seq_len(n)), rep(seq_len(K), length.out = n))
  mean_SSR <- matrix(NA_real_, nrow = length(lambda_a_seq), ncol = length(lambda_b_seq))
  dimnames(mean_SSR) <- list(
    paste0("la=", lambda_a_seq),
    paste0("lb=", lambda_b_seq)
  )
  
  for (ia in seq_along(lambda_a_seq)) {
    for (ib in seq_along(lambda_b_seq)) {
      lam_a <- lambda_a_seq[ia]
      lam_b <- lambda_b_seq[ib]
      cv_SSR <- numeric(K)
      
      for (k in seq_len(K)) {
        test_idx <- folds[[k]]
        train_idx <- setdiff(seq_len(n), test_idx)
        
        X_train <- X[train_idx, , drop = FALSE]
        M_train <- M[train_idx, , drop = FALSE]
        Y_train <- Y[train_idx]
        X_test <- X[test_idx, , drop = FALSE]
        M_test <- M[test_idx, , drop = FALSE]
        Y_test <- Y[test_idx]
        
        model <- optimize_weights(
          X_train, M_train, Y_train,
          lambda_n = lambda_n,
          lambda_a = lam_a,
          lambda_b = lam_b,
          max_iter = 50,
          tol = 1e-3
        )
        
        a_hat <- model$a
        b_hat <- model$b
        x_test <- as.vector(X_test %*% a_hat)
        z_test <- as.vector(M_test %*% b_hat)
        test_fit <- lm(Y_test ~ x_test + z_test)
        cv_SSR[k] <- sum(test_fit$residuals^2)
      }
      
      mean_SSR[ia, ib] <- mean(cv_SSR)
    }
  }
  
  min_idx <- which(mean_SSR == min(mean_SSR, na.rm = TRUE), arr.ind = TRUE)[1, ]
  list(
    best_lambda_a = lambda_a_seq[min_idx[1]],
    best_lambda_b = lambda_b_seq[min_idx[2]],
    cv_error_grid = mean_SSR
  )
}

fit_memm_on_dataset <- function(sim_data,
                                lambda_n = 1,
                                lambda_a_seq = c(0, 0.1, 0.2, 0.5),
                                lambda_b_seq = c(0, 0.1, 0.2, 0.5),
                                K = 5,
                                final_n_restarts = 1,
                                restart_seed = NULL) {
  X <- sim_data$X
  M <- sim_data$M
  Y <- sim_data$Y
  
  cv_result <- cv_select_lambda(
    X = X,
    M = M,
    Y = Y,
    lambda_n = lambda_n,
    lambda_a_seq = lambda_a_seq,
    lambda_b_seq = lambda_b_seq,
    K = K
  )
  
  model <- fit_with_restarts(
    X = X,
    M = M,
    Y = Y,
    lambda_n = lambda_n,
    lambda_a = cv_result$best_lambda_a,
    lambda_b = cv_result$best_lambda_b,
    n_restarts = final_n_restarts,
    max_iter = 50,
    tol = 1e-3,
    restart_seed = restart_seed
  )
  
  list(
    a = model$a,
    b = model$b,
    best_lambda_a = cv_result$best_lambda_a,
    best_lambda_b = cv_result$best_lambda_b
  )
}

build_all_method_matrix_list <- function(
    sim_data,
    methods = c("MEMM", "FBAS", "SIS_MCP", "DM1",
                "DM2", "SPCMA", "Pathway", "TS"),
    threshold = 1e-2,
    lambda_fbas = 0.01,
    rho_fbas = 1,
    lambda_n = 1,
    lambda_a_seq = c(0, 0.1, 0.2, 0.5),
    lambda_b_seq = c(0, 0.1, 0.2, 0.5),
    K = 5,
    final_n_restarts = 1,
    restart_seed = NULL) {
  out <- list()
  comparison_methods <- setdiff(methods, "MEMM")
  
  if (length(comparison_methods) > 0) {
    out <- build_method_b_list(
      sim_data = sim_data,
      threshold = threshold,
      lambda_fbas = lambda_fbas,
      rho_fbas = rho_fbas,
      methods = comparison_methods
    )
  }
  
  if ("MEMM" %in% methods) {
    memm_fit <- fit_memm_on_dataset(
      sim_data = sim_data,
      lambda_n = lambda_n,
      lambda_a_seq = lambda_a_seq,
      lambda_b_seq = lambda_b_seq,
      K = K,
      final_n_restarts = final_n_restarts,
      restart_seed = restart_seed
    )
    out <- c(list(MEMM = tcrossprod(memm_fit$a, memm_fit$b)), out)
  }
  
  out
}

build_table12_design <- function(
    pathways = c("complete", "partial"),
    size_grid = data.frame(
      m = c(20L, 50L),
      q = c(20L, 50L),
      stringsAsFactors = FALSE
    ),
    rhoX_values = c(0, 0.3),
    rhoM_values = c(0.1, 0.2),
    n1 = 200L) {
  if (!all(c("m", "q") %in% names(size_grid))) {
    stop("size_grid must contain columns m and q.")
  }
  
  base_design <- merge(
    expand.grid(
      pathway = pathways,
      size_idx = seq_len(nrow(size_grid)),
      rhoX = rhoX_values,
      rhoM = rhoM_values,
      stringsAsFactors = FALSE
    ),
    cbind(size_idx = seq_len(nrow(size_grid)), size_grid),
    by = "size_idx",
    sort = FALSE
  )
  base_design$n1 <- n1
  base_design$q_a <- pmax(1L, round(base_design$q * 0.25))
  base_design$r <- pmax(1L, round(base_design$m * 0.25))
  base_design <- base_design[order(base_design$pathway, base_design$m,
                                   base_design$q, base_design$rhoX,
                                   base_design$rhoM), , drop = FALSE]
  base_design$scenario_id <- ave(
    base_design$pathway,
    base_design$pathway,
    FUN = function(x) seq_along(x)
  )
  rownames(base_design) <- NULL
  base_design
}

run_method_table_replications <- function(
    nrep = 100,
    n1 = 200,
    m = 20,
    q = 20,
    q_a = max(1L, round(q * 0.25)),
    r = max(1L, round(m * 0.25)),
    rhoX = 0,
    rhoM = 0.1,
    sigma1 = 1,
    sigma2 = 1,
    pathway = c("complete", "partial", "none"),
    x_dist = c("normal", "t"),
    m_noise_dist = c("normal", "t"),
    df_x = 6,
    df_m = 6,
    c_signal = NULL,
    methods = c("MEMM", "FBAS", "SIS_MCP", "DM1",
                "DM2", "SPCMA", "Pathway", "TS"),
    threshold = 1e-2,
    lambda_fbas = 0.01,
    rho_fbas = 1,
    lambda_n = 1,
    lambda_a_seq = c(0, 0.1, 0.2, 0.5),
    lambda_b_seq = c(0, 0.1, 0.2, 0.5),
    K = 5,
    final_n_restarts = 1,
    restart_seed = NULL,
    proxy_aggregation = c("l2", "l1"),
    truth_mode = c("legacy_mediator_only", "outer_active"),
    seed = NULL,
    verbose = TRUE) {
  pathway <- match.arg(pathway)
  x_dist <- match.arg(x_dist)
  m_noise_dist <- match.arg(m_noise_dist)
  proxy_aggregation <- match.arg(proxy_aggregation)
  truth_mode <- match.arg(truth_mode)
  
  out <- vector("list", nrep)
  for (rep_idx in seq_len(nrep)) {
    if (!is.null(seed)) {
      set.seed(seed + rep_idx - 1L)
    }
    if (verbose && (rep_idx == 1L || rep_idx %% 10L == 0L || rep_idx == nrep)) {
      cat(sprintf("  replication %d/%d\n", rep_idx, nrep))
    }
    
    sim_data <- simulate_data(
      n1 = n1,
      m = m,
      q = q,
      q_a = q_a,
      r = r,
      rhoX = rhoX,
      rhoM = rhoM,
      sigma1 = sigma1,
      sigma2 = sigma2,
      pathway = pathway,
      x_dist = x_dist,
      m_noise_dist = m_noise_dist,
      df_x = df_x,
      df_m = df_m,
      c_signal = c_signal
    )
    
    rep_restart_seed <- if (is.null(restart_seed)) NULL else restart_seed + rep_idx - 1L
    method_matrix_list <- build_all_method_matrix_list(
      sim_data = sim_data,
      methods = methods,
      threshold = threshold,
      lambda_fbas = lambda_fbas,
      rho_fbas = rho_fbas,
      lambda_n = lambda_n,
      lambda_a_seq = lambda_a_seq,
      lambda_b_seq = lambda_b_seq,
      K = K,
      final_n_restarts = final_n_restarts,
      restart_seed = rep_restart_seed
    )
    
    res_i <- evaluate_method_matrices(
      method_matrix_list = method_matrix_list,
      true_a = sim_data$true_a,
      true_b = sim_data$true_b,
      true_active_mediators = sim_data$active_mediators,
      threshold = threshold,
      proxy_aggregation = proxy_aggregation,
      truth_mode = truth_mode
    )
    res_i$replication <- rep_idx
    out[[rep_idx]] <- res_i
  }
  
  do.call(rbind, out)
}

run_table12_summaries <- function(
    nrep = 100,
    pathways = c("complete", "partial"),
    size_grid = data.frame(
      m = c(20L, 50L),
      q = c(20L, 50L),
      stringsAsFactors = FALSE
    ),
    rhoX_values = c(0, 0.3),
    rhoM_values = c(0.1, 0.2),
    n1 = 200,
    sigma1 = 1,
    sigma2 = 1,
    x_dist = "normal",
    m_noise_dist = "normal",
    df_x = 6,
    df_m = 6,
    c_signal = NULL,
    methods = c("MEMM", "FBAS", "SIS_MCP", "DM1",
                "DM2", "SPCMA", "Pathway", "TS"),
    threshold = 1e-2,
    lambda_fbas = 0.01,
    rho_fbas = 1,
    lambda_n = 1,
    lambda_a_seq = c(0, 0.1, 0.2, 0.5),
    lambda_b_seq = c(0, 0.1, 0.2, 0.5),
    K = 5,
    final_n_restarts = 1,
    restart_seed = NULL,
    proxy_aggregation = "l2",
    truth_mode = "legacy_mediator_only",
    seed = NULL,
    verbose = TRUE) {
  design <- build_table12_design(
    pathways = pathways,
    size_grid = size_grid,
    rhoX_values = rhoX_values,
    rhoM_values = rhoM_values,
    n1 = n1
  )
  
  scenario_results <- vector("list", nrow(design))
  for (i in seq_len(nrow(design))) {
    scenario_i <- design[i, , drop = FALSE]
    if (verbose) {
      cat(sprintf(
        "\nScenario %d/%d: pathway=%s, m=%d, q=%d, rhoX=%.2f, rhoM=%.2f\n",
        i, nrow(design), scenario_i$pathway, scenario_i$m, scenario_i$q,
        scenario_i$rhoX, scenario_i$rhoM
      ))
    }
    
    scenario_seed <- if (is.null(seed)) NULL else seed + (i - 1L) * nrep
    scenario_restart_seed <- if (is.null(restart_seed)) NULL else restart_seed + (i - 1L) * nrep
    res_i <- run_method_table_replications(
      nrep = nrep,
      n1 = scenario_i$n1,
      m = scenario_i$m,
      q = scenario_i$q,
      q_a = scenario_i$q_a,
      r = scenario_i$r,
      rhoX = scenario_i$rhoX,
      rhoM = scenario_i$rhoM,
      sigma1 = sigma1,
      sigma2 = sigma2,
      pathway = scenario_i$pathway,
      x_dist = x_dist,
      m_noise_dist = m_noise_dist,
      df_x = df_x,
      df_m = df_m,
      c_signal = c_signal,
      methods = methods,
      threshold = threshold,
      lambda_fbas = lambda_fbas,
      rho_fbas = rho_fbas,
      lambda_n = lambda_n,
      lambda_a_seq = lambda_a_seq,
      lambda_b_seq = lambda_b_seq,
      K = K,
      final_n_restarts = final_n_restarts,
      restart_seed = scenario_restart_seed,
      proxy_aggregation = proxy_aggregation,
      truth_mode = truth_mode,
      seed = scenario_seed,
      verbose = FALSE
    )
    res_i$pathway <- scenario_i$pathway
    res_i$scenario_id <- scenario_i$scenario_id
    res_i$m <- scenario_i$m
    res_i$q <- scenario_i$q
    res_i$rhoX <- scenario_i$rhoX
    res_i$rhoM <- scenario_i$rhoM
    scenario_results[[i]] <- res_i
  }
  
  combined <- do.call(rbind, scenario_results)
  summary_df <- summarize_method_table(
    combined,
    group_cols = c("pathway", "scenario_id", "m", "q", "rhoX", "rhoM", "method"),
    include_sd = FALSE
  )
  summary_df <- format_table12_summary(summary_df)
  
  list(
    design = design,
    raw_results = combined,
    summary = summary_df,
    table1_complete = subset(summary_df, pathway == "complete"),
    table2_partial = subset(summary_df, pathway == "partial")
  )
}

###############################################################################
# 5) Original simulation wrapper, plus lightweight summary columns
###############################################################################

run_simulation_with_cv <- function(n_runs = 10,
                                   n1 = 100, m = 20, q = 10, q_a = 3, r = 3,
                                   rhoX = 0.3, rhoM = 0.3,
                                   sigma1 = 1, sigma2 = 1,
                                   pathway = c("partial", "complete", "none"),
                                   x_dist = c("normal", "t"),
                                   m_noise_dist = c("normal", "t"),
                                   df_x = 6,
                                   df_m = 6,
                                   c_signal = NULL,
                                   lambda_n = 1,
                                   lambda_a_seq = c(0, 0.1, 0.2, 0.5),
                                   lambda_b_seq = c(0, 0.1, 0.2, 0.5),
                                   K = 5,
                                   threshold = 1e-3,
                                   final_n_restarts = 1,
                                   restart_seed = NULL,
                                   seed = NULL,
                                   verbose = TRUE) {
  pathway <- match.arg(pathway)
  x_dist <- match.arg(x_dist)
  m_noise_dist <- match.arg(m_noise_dist)
  
  true_MP <- true_mp_from_pathway(pathway)
  results <- data.frame(
    MP = numeric(0),
    true_MP = numeric(0),
    AbsBias_MP = numeric(0),
    Accuracy = numeric(0),
    Precision = numeric(0),
    Recall = numeric(0),
    F1 = numeric(0),
    L2_error_a = numeric(0),
    L2_error_b = numeric(0),
    Cosine_a = numeric(0),
    Cosine_b = numeric(0),
    Cosine_avg = numeric(0)
  )
  
  for (run in seq_len(n_runs)) {
    if (!is.null(seed)) {
      set.seed(seed + run - 1L)
    }
    
    data <- simulate_data(
      n1 = n1, m = m, q = q, q_a = q_a, r = r,
      rhoX = rhoX, rhoM = rhoM,
      sigma1 = sigma1, sigma2 = sigma2,
      pathway = pathway,
      x_dist = x_dist,
      m_noise_dist = m_noise_dist,
      df_x = df_x,
      df_m = df_m,
      c_signal = c_signal
    )
    
    X <- data$X
    M <- data$M
    Y <- data$Y
    true_active_med <- data$active_mediators
    
    cv_result <- cv_select_lambda(
      X, M, Y, lambda_n,
      lambda_a_seq, lambda_b_seq,
      K = K
    )
    
    run_restart_seed <- if (is.null(restart_seed)) NULL else restart_seed + run - 1L
    
    model <- fit_with_restarts(
      X, M, Y,
      lambda_n = lambda_n,
      lambda_a = cv_result$best_lambda_a,
      lambda_b = cv_result$best_lambda_b,
      n_restarts = final_n_restarts,
      max_iter = 50,
      tol = 1e-3,
      restart_seed = run_restart_seed
    )
    a_hat <- model$a
    b_hat <- model$b
    
    ## Summarize final model
    x_agg <- as.vector(X %*% a_hat)
    z_agg <- as.vector(M %*% b_hat)
    fit <- lm(Y ~ x_agg + z_agg)
    coef_vals <- coef(fit)
    gamma_hat <- as.numeric(coef_vals["x_agg"])
    eta_hat <- as.numeric(coef_vals["z_agg"])
    alpha_hat <- as.numeric(crossprod(x_agg, z_agg))
    denom <- alpha_hat * eta_hat + gamma_hat
    MP_hat <- if (abs(denom) < 1e-8) 0 else (alpha_hat * eta_hat / denom)
    MP_hat <- max(0, min(1, MP_hat))
    
    direction_a <- compute_direction_metrics(a_hat, data$true_a)
    direction_b <- compute_direction_metrics(b_hat, data$true_b)
    
    ## Identify active mediators only, as in the original first draft
    pred_active_med <- which(abs(b_hat) > threshold)
    true_active_set <- if (length(true_active_med) > 0) true_active_med else integer(0)
    
    TP <- length(intersect(pred_active_med, true_active_set))
    FP <- length(setdiff(pred_active_med, true_active_set))
    FN <- length(setdiff(true_active_set, pred_active_med))
    TN <- q - (TP + FP + FN)
    
    if ((TP + FP + FN) == 0) {
      accuracy <- 1
      precision <- 1
      recall <- 1
      F1 <- 1
    } else {
      accuracy <- (TP + TN) / q
      precision <- if ((TP + FP) > 0) TP / (TP + FP) else 0
      recall <- if ((TP + FN) > 0) TP / (TP + FN) else 0
      if ((precision + recall) == 0) {
        F1 <- 0
      } else {
        F1 <- 2 * precision * recall / (precision + recall)
      }
    }
    
    results <- rbind(
      results,
      data.frame(
        MP = MP_hat,
        true_MP = true_MP,
        AbsBias_MP = abs(MP_hat - true_MP),
        Accuracy = accuracy,
        Precision = precision,
        Recall = recall,
        F1 = F1,
        L2_error_a = direction_a$l2_error,
        L2_error_b = direction_b$l2_error,
        Cosine_a = direction_a$cosine,
        Cosine_b = direction_b$cosine,
        Cosine_avg = compute_average_cosine(direction_a$cosine, direction_b$cosine)
      )
    )
  }
  
  if (verbose) {
    avg_metrics <- colMeans(
      results[, c("Accuracy", "Precision", "Recall", "F1",
                  "L2_error_a", "L2_error_b", "Cosine_a", "Cosine_b", "Cosine_avg")],
      na.rm = TRUE
    )
    avg_MP <- mean(results$MP, na.rm = TRUE)
    avg_abs_bias <- mean(results$AbsBias_MP, na.rm = TRUE)
    
    cat(sprintf("Avg MP: %.3f\n", avg_MP))
    cat(sprintf("Avg |MP bias|: %.3f\n", avg_abs_bias))
    cat(sprintf("Avg Accuracy: %.3f\n", avg_metrics["Accuracy"]))
    cat(sprintf("Avg Precision: %.3f\n", avg_metrics["Precision"]))
    cat(sprintf("Avg Recall: %.3f\n", avg_metrics["Recall"]))
    cat(sprintf("Avg F1-score: %.3f\n", avg_metrics["F1"]))
    cat(sprintf("Avg L2 error (a): %.3f\n", avg_metrics["L2_error_a"]))
    cat(sprintf("Avg L2 error (b): %.3f\n", avg_metrics["L2_error_b"]))
    cat(sprintf("Avg cosine (a): %.3f\n", avg_metrics["Cosine_a"]))
    cat(sprintf("Avg cosine (b): %.3f\n", avg_metrics["Cosine_b"]))
    cat(sprintf("Avg cosine (avg): %.3f\n", avg_metrics["Cosine_avg"]))
  }
  
  results
}

###############################################################################
# 6) Summaries and scenario grids
###############################################################################

summarize_simulation_results <- function(results_df) {
  results_df <- add_average_cosine_columns(results_df)
  metric_cols <- c(
    "MP", "true_MP", "AbsBias_MP",
    "Accuracy", "Precision", "Recall", "F1",
    "L2_error_a", "L2_error_b", "Cosine_a", "Cosine_b", "Cosine_avg"
  )
  metric_cols <- intersect(metric_cols, names(results_df))
  
  mean_vals <- stats::setNames(
    vapply(metric_cols, function(col) mean(results_df[[col]], na.rm = TRUE), numeric(1)),
    paste0(metric_cols, "_mean")
  )
  sd_vals <- stats::setNames(
    vapply(metric_cols, function(col) stats::sd(results_df[[col]], na.rm = TRUE), numeric(1)),
    paste0(metric_cols, "_sd")
  )
  
  as.data.frame(as.list(c(mean_vals, sd_vals)), check.names = FALSE)
}

run_scenario_grid <- function(design_grid, common_args,
                              seed = NULL, verbose = TRUE) {
  if (!is.data.frame(design_grid) || nrow(design_grid) == 0) {
    stop("design_grid must be a non-empty data.frame.")
  }
  if (!is.list(common_args)) {
    stop("common_args must be a list.")
  }
  
  out <- vector("list", nrow(design_grid))
  valid_names <- names(formals(run_simulation_with_cv))
  
  for (i in seq_len(nrow(design_grid))) {
    row_args <- as.list(design_grid[i, , drop = FALSE])
    row_args <- row_args[intersect(names(row_args), valid_names)]
    if (length(row_args) > 0) {
      keep_row_args <- !vapply(
        row_args,
        function(x) length(x) == 1 && is.atomic(x) && is.na(x),
        logical(1)
      )
      row_args <- row_args[keep_row_args]
    }
    call_args <- utils::modifyList(common_args, row_args)
    call_args$verbose <- FALSE
    if (!is.null(seed)) {
      call_args$seed <- seed + i - 1L
    }
    
    fit_i <- do.call(run_simulation_with_cv, call_args)
    sum_i <- summarize_simulation_results(fit_i)
    out[[i]] <- cbind(design_grid[i, , drop = FALSE], sum_i, row.names = NULL)
    
    if (verbose) {
      tag_cols <- intersect(c("scenario", "pathway", "rhoX", "rhoM", "c_signal", "df_x", "df_m"),
                            names(design_grid))
      tag_msg <- paste(
        vapply(tag_cols,
               function(col) sprintf("%s=%s", col, as.character(design_grid[i, col])),
               character(1)),
        collapse = ", "
      )
      cat(sprintf("Completed scenario %d/%d: %s\n", i, nrow(design_grid), tag_msg))
    }
  }
  
  do.call(rbind, out)
}

    

