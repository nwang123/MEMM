MEMM
===
This R package provides a simulation and estimation framework for high-dimensional multivariate mediation analysis, integrating exposure, mediator, and outcome data through a penalized joint model solved by ADMM (Alternating Direction Method of Multipliers).
It supports realistic simulation of correlated omics-style data, cross-validated tuning of regularization parameters, and evaluation of mediation performance metrics such as accuracy, precision, recall, F1-score, and estimated mediation proportion (MP).

![Alt text](./Fig1.png)

## R Code Overview

The implementation consists of four main functions that together form a complete pipeline for simulation, model estimation, and performance evaluation.
	1.	simulate_data()
Generates synthetic exposure–mediator–outcome data under user-specified correlation, noise, and mediation structures.
It allows for complete, partial, or no mediation pathways by modifying the true parameters (\alpha, \eta, \gamma).
The function returns standardized matrices X, M, Y and true coefficients for downstream benchmarking.
	2.	optimize_weights()
Implements the ADMM-based algorithm to estimate the projection directions a and b under joint LASSO penalties \lambda_a and \lambda_b, with additional penalty \lambda_n for the correlation constraint.
Iteratively updates primal and dual variables until convergence, returning the estimated a, b.
	3.	cv_select_lambda()
Performs K-fold cross-validation over a grid of (\lambda_a, \lambda_b) to minimize the average sum of squared residuals (SSR).
Outputs the optimal tuning parameters and the full error surface for diagnostic visualization.
	4.	run_simulation_with_cv()
Serves as the master wrapper that integrates data generation, cross-validation, ADMM fitting, and performance evaluation (Accuracy, Precision, Recall, F1, and Mediation Proportion).
Designed for repeated Monte Carlo runs to summarize mean performance across simulation replicates.

Installation
===
To install the MEMM package, you will first need to install devtools package and then execute the following code:
```
#install.packages("devtools")
library(devtools)
install_github("nwang123/MEMM")
```
Usage
===========
The following help page will also provide quick references for TIPS package and the example command lines:
```
library(MEMM)
```

Output 
===========
The main output of our tool is a list of pathways and genes with p-values from likelihood ratio test for their association with the phenotype of interest. These findings could offer more reliable and interpretable results for TWAS analyses. 

Development
===========
This R package is developed by Neng Wang.
