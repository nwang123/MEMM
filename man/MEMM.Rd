\name{memm}
\alias{memm}
\title{Estimation framework for high-dimensional multivariate mediation analysis, integrating exposure, mediator, and outcome data through a penalized joint model solved by ADMM (Alternating Direction Method of Multipliers)}
\description{
The MEMM package provides a comprehensive framework for high-dimensional multivariate mediation analysis. It implements a penalized joint model that integrates exposure, mediator, and outcome data to estimate mediation effects in high-dimensional settings. The model is solved using the Alternating Direction Method of Multipliers (ADMM) algorithm with cross-validation for parameter selection. The package includes functions for data simulation, model estimation, cross-validation, and comprehensive performance evaluation including mediation proportion estimation and variable selection metrics.
}
\usage{
simulate_data(n1, m, q, q_a, r, rhoX = 0, rhoM = 0, sigma1 = 1, sigma2 = 1, 
              pathway = c("complete","partial","none"))
              
optimize_weights(X, M, Y, lambda_n, lambda_a, lambda_b, max_iter = 50, tol = 1e-4)

cv_select_lambda(X, M, Y, lambda_n, lambda_a_seq, lambda_b_seq, K = 5)

run_simulation_with_cv(n_runs = 10, n1 = 100, m = 20, q = 10, q_a = 3, r = 3,
                       rhoX = 0.3, rhoM = 0.3, sigma1 = 1, sigma2 = 1,
                       pathway = c("partial","complete","none"), lambda_n = 1,
                       lambda_a_seq = c(0, 0.1, 0.2, 0.5),
                       lambda_b_seq = c(0, 0.1, 0.2, 0.5), K = 5, threshold = 1e-3)
}
\arguments{
  \item{n1}{Number of observations}
  \item{m}{Number of exposures}
  \item{q}{Number of mediators}
  \item{q_a}{Number of active mediators (with non-zero effects)}
  \item{r}{Number of active exposures (with non-zero effects)}
  \item{rhoX}{Correlation parameter for exposure matrix (default: 0)}
  \item{rhoM}{Correlation parameter for mediator matrix (default: 0)}
  \item{sigma1}{Noise factor for mediators (default: 1)}
  \item{sigma2}{Noise factor for outcome (default: 1)}
  \item{pathway}{Type of mediation pathway: "complete" (pure mediation), "partial" (partial mediation), or "none" (no mediation)}
  \item{X}{Exposure matrix (n1 x m)}
  \item{M}{Mediator matrix (n1 x q)}
  \item{Y}{Outcome vector (length n1)}
  \item{lambda_n}{Tuning parameter for non-convex penalty}
  \item{lambda_a}{Tuning parameter for L1 penalty on exposure weights}
  \item{lambda_b}{Tuning parameter for L1 penalty on mediator weights}
  \item{max_iter}{Maximum number of ADMM iterations (default: 50)}
  \item{tol}{Convergence tolerance (default: 1e-4)}
  \item{lambda_a_seq}{Sequence of lambda_a values for cross-validation}
  \item{lambda_b_seq}{Sequence of lambda_b values for cross-validation}
  \item{K}{Number of folds for cross-validation (default: 5)}
  \item{n_runs}{Number of simulation runs (default: 10)}
  \item{threshold}{Threshold for identifying active mediators (default: 1e-3)}
}
\details{
The MEMM framework addresses the challenges of high-dimensional mediation analysis by:

1. Using a joint modeling approach that simultaneously estimates exposure and mediator effects
2. Implementing L1 regularization for sparse solutions in high-dimensional settings
3. Employing ADMM for efficient optimization with convergence guarantees
4. Incorporating cross-validation for automatic parameter selection
5. Providing comprehensive simulation tools for method evaluation

The \code{simulate_data} function generates realistic multivariate data with specified mediation structures. The \code{optimize_weights} function implements the core ADMM algorithm. The \code{cv_select_lambda} function performs cross-validation for parameter tuning, and \code{run_simulation_with_cv} provides a complete simulation workflow with performance evaluation.
}
\value{

}
\examples{
}
\seealso{
Useful links:
\itemize{
  \item \pkg{glmnet} for alternative high-dimensional regression approaches
  \item \pkg{mediation} for traditional mediation analysis
  \item \pkg{HIMA} for high-dimensional mediation analysis
}
}
\author{
Neng Wang
}
