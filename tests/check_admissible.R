# Run from the repository root: Rscript tests/check_admissible.R
source("R/MEMM.R")
set.seed(20261008)
X <- scale(matrix(rnorm(80*5),80,5),scale=FALSE)
M <- scale(matrix(rnorm(80*6),80,6)+X[,1],scale=FALSE)
Y <- as.vector(scale(0.5*X[,1]+0.8*M[,1]+rnorm(80),scale=FALSE))
initial <- memm_initialize(X,M,Y,NULL,NULL,1e-6,1e-6)
a <- initial$a; b <- initial$b
x <- as.vector(X %*% a); m <- as.vector(M %*% b)
yx <- lm(Y ~ x - 1); mx <- lm(m ~ x - 1); yxm <- lm(Y ~ x + m - 1)
alpha <- unname(coef(mx)); tau <- unname(coef(yx)); eta <- unname(coef(yxm)[2])
reference <- (sum(resid(yx)^2)+sum(resid(mx)^2)+
                sum(resid(yxm)^2)/(1-alpha^2))/(2*nrow(X)) -
  0.1*alpha*eta/tau + 0.02*sum(abs(a)) + 0.03*sum(abs(b))
stopifnot(abs(reference-memm_objective(X,M,Y,a,b,0.1,0.02,0.03)) < 1e-10)
cat("PASS: Eq. (4) agrees with independent OLS residual calculations.\n")

# Check analytic derivatives against finite differences along the NORMALIZED
# loading path. This detects missing profiling derivatives and SSR scaling.
gradient <- memm_smooth_gradient(X,M,Y,a,b,0.1)
for (block in c("a","b")) {
  weight <- if (block=="a") a else b
  design <- if (block=="a") X else M
  tangent_normal <- as.vector(crossprod(design,design %*% weight))
  for (j in 1:5) {
    v <- rnorm(length(weight))
    v <- v-tangent_normal*sum(tangent_normal*v)/sum(tangent_normal^2)
    v <- v/sqrt(sum(v^2))
    evaluate <- function(h) {
      w <- memm_normalize(weight+h*v,design)
      memm_objective(X,M,Y,if(block=="a")w else a,if(block=="b")w else b,0.1,0,0)
    }
    h <- 1e-6
    fd <- (evaluate(h)-evaluate(-h))/(2*h)
    analytic <- sum(gradient[[block]]*v)
    stopifnot(abs(fd-analytic) < 1e-5*max(1,abs(fd),abs(analytic)))
  }
}
cat("PASS: both profiled loading gradients agree with finite differences.\n")

fit <- optimize_weights(X,M,Y,0.1,0.02,0.03,max_iter=30,inner_max_iter=6)
h <- fit$history
stopifnot(all(h$tau_margin>=0),all(h$alpha_margin>=0),
          max(h$norm_error_a,h$norm_error_b)<1e-8,
          fit$diagnostics$accepted_feasibility_violations==0,
          !fit$diagnostics$theorem_conditions_verified)
for (i in which(h$stage %in% c("after_a","after_b"))) {
  stopifnot(h$augmented_lagrangian[i] <= h$augmented_lagrangian[i-1]+1e-10)
}
stopifnot(max(abs(fit$u_a))<=0.02+1e-10,max(abs(fit$u_b))<=0.03+1e-10)
stopifnot(isTRUE(all.equal(fit$objective,memm_objective(X,M,Y,fit$a,fit$b,0.1,0.02,0.03))))
cat("PASS: intermediate/complete feasibility, block decrease and scaled dual bounds.\n")

# Force a collinear initial pair and a separation bound to exercise feasibility
# repair and rejected trial steps, rather than only testing interior iterates.
Z <- cbind(X[,1],X[,2],X[,3])
boundary <- optimize_weights(X,Z,Y,0.2,0.02,0.03,max_iter=8,inner_max_iter=5,
  init_a=c(1,0,0,0,0),init_b=c(1,0,0),delta=0.9,step_size=100)
stopifnot(boundary$diagnostics$initialization_repaired["b"],
          boundary$diagnostics$rejected_feasibility_trials>0,
          all(boundary$history$alpha_margin>=0),all(boundary$history$tau_margin>=0))
bad <- try(optimize_weights(X,M,Y,0.1,0.02,0.03,r0=1e10),silent=TRUE)
stopifnot(inherits(bad,"try-error"))
# Rank-deficient loading spaces still require feasible aggregates.
rankdef <- optimize_weights(cbind(X,X),cbind(M,M),Y,0.1,0.02,0.03,
                           max_iter=4,inner_max_iter=3)
stopifnot(all(rankdef$history$tau_margin>=0),all(rankdef$history$alpha_margin>=0))
cat("PASS: inadmissible initialization, infeasible r0, boundary rejection and rank deficiency.\n")

restart <- fit_with_restarts(X,M,Y,0.1,0.02,0.03,n_restarts=2,max_iter=5,
  restart_seed=17,optimizer_control=list(delta=0.2,inner_max_iter=3))
stopifnot(restart$best_score==min(restart$restart_scores,na.rm=TRUE),restart$control$delta==0.2)
cv <- cv_select_lambda(X,M,Y,0.1,0.02,0.03,K=2,
                       optimizer_control=list(inner_max_iter=2,delta=0.2))
stopifnot(all(is.finite(cv$cv_error_grid)))
cat("PASS: restart selection uses Eq. (4); controls reach cross-validation.\n")
cat("These finite checks do not establish Theorem 1's sequence/geometric assumptions.\n")
