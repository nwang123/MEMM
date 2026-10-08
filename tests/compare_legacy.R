# A diagnostic comparison, not a rerun of the manuscript's simulation study.
source("R/MEMM.R")
path <- commandArgs(trailingOnly=TRUE)[1]
legacy <- new.env(parent=globalenv())
for (expr in parse(path)) {
  if (is.call(expr) && identical(expr[[1]],as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name("function"))) {
    eval(expr,legacy)
  }
}
set.seed(20261008)
X <- scale(matrix(rnorm(100*5),100,5),scale=FALSE)
M <- scale(matrix(rnorm(100*6),100,6)+X[,1],scale=FALSE)
Y <- as.vector(scale(0.5*X[,1]+0.8*M[,1]+rnorm(100),scale=FALSE))
old <- legacy$optimize_weights(X,M,Y,0.1,0.02,0.03,max_iter=50,tol=1e-4)
new <- optimize_weights(X,M,Y,0.1,0.02,0.03,max_iter=50,tol=1e-4)
summarize <- function(fit,label) {
  p <- memm_profile(X,M,Y,fit$a,fit$b)
  f <- memm_feasibility(X,M,Y,fit$a,fit$b)
  data.frame(version=label,MP=p$MP,tau=p$tau,alpha=p$alpha,
             feasible=f$feasible,tau_margin=f$tau_margin,alpha_margin=f$alpha_margin,
             eq4=memm_objective(X,M,Y,fit$a,fit$b,0.1,0.02,0.03))
}
print(rbind(summarize(old,"legacy"),summarize(new,"revised")),row.names=FALSE,digits=8)
print(new$diagnostics)
q <- new$history$decrease_ratio[new$history$stage=="complete"]
cat("Iterations with negative observed primal decrease ratio:",sum(q<0,na.rm=TRUE),"\n")
cat("Minimum observed primal decrease ratio:",min(q,na.rm=TRUE),"\n")
complete_objectives <- new$history$objective[new$history$stage %in% c("initial","complete")]
cat("Maximum complete-iteration primal objective increase:",max(diff(complete_objectives)),"\n")
cat("No equivalence or negligible-impact claim follows from this example.\n")
