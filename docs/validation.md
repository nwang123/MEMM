# Validation of the admissible Equation (4) revision

Status: draft correction, **not ready for a claim of unchanged manuscript results**.

Code tested: commit `1b58d13b250b4904ab4fa850c3a4202f25cec5b2`.
GitHub Actions: https://github.com/nwang123/MEMM/actions/runs/37785200901
R 4.3.3, Ubuntu; all automated code checks passed.

## Numerical checks

The objective agrees with independent OLS residual calculations. Both loading
gradients agree with central finite differences along normalized loading paths.
Tests also passed for intermediate feasibility, block augmented-objective
decrease, scaled dual bounds, invalid initialization, impossible r0, rejection
near the separation boundary, rank-deficient designs, and propagation of
controls through cross-validation and restarts. These are implementation checks,
not verification of Theorem 1's assumptions.

## Fixed synthetic comparison

`tests/compare_legacy.R` uses the same centered data, seed, tuning values,
maximum iterations and outer tolerance for both code versions. It does not
reproduce the manuscript's simulation design.

| Version | Signed MP | Estimated total effect | Equation (4) |
|---|---:|---:|---:|
| Legacy optimizer | 1.08348887 | 8.9398853 | 4.9717075 |
| Revised optimizer | 0.16548395 | 9.8195488 | 1.2862040 |

Both final pairs were feasible. The revised run rejected 8 infeasible trials
and accepted no infeasible loading pairs. Its minimum accepted feasibility
margins were 0.2012531 for tau-r0 and 0.235699 for 1-alpha^2-delta.
The final primal/dual infinity-norm residuals were approximately 3.47e-18 and
8.43e-05. One complete iteration increased the primal objective by
2.540162e-07; the minimum observed decrease ratio was -0.003850781.
Thus the observed path cannot be reported as verification of (B2).

This comparison changes the objective-gradient implementation and numerical
block solver as well as enforcing D. It does not isolate the effect of the
feasibility guard alone. No negligible-impact conclusion is justified.

## Scope and next step

This workflow reports only synthetic-data diagnostics. Reassess the actual
manuscript analyses before claiming that the implementation changes preserve
its findings. The numerical choices of r0 and delta need to be justified for
the analysis scale: a small positive r0 prevents division by zero but does
not mathematically rule out a very large signed MP. Numerical stopping and
feasibility checks do not verify Theorem 1's sequence or geometric assumptions.
