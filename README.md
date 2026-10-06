# tppcm

<!-- badges: start -->
[![R-CMD-check](https://github.com/jiewenTsai/tppcm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jiewenTsai/tppcm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Boundary-level diagnostics for the two-parameter partial credit model
(TPPCM; Yu, 1991), built on top of **TAM**.

In the TPPCM every step *l* of item *i* has its own discrimination
*a<sub>il</sub>*:

P(X<sub>i</sub> = k | θ) ∝ exp( Σ<sub>l ≤ k</sub> (a<sub>il</sub> θ − d<sub>il</sub>) ),  k = 0, …, K−1.

The generalized partial credit model (GPCM; Muraki, 1992) is the special case
with equal steps within items, and the partial credit model (PCM; Masters,
1982) has one common discrimination.

Estimation stays in TAM (and in sirt for the product form). tppcm takes the
fitted object and adds what TAM does not report at the level of category
boundaries: step discriminations (CBD), locations (CBL) and intercepts (IBD)
with standard errors, tests between models, fit statistics, step information
and tests of differential step functioning (DSF) across groups.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jiewenTsai/tppcm", build_vignettes = TRUE)
```

## Features

* **Step designs for TAM.** `tppcm(dat, design = )` builds the design array
  for `tam.mml.3pl()`. The design formula reads the step discriminations as
  an item x step layout on the log scale: `~ item * step` (saturated TPPCM,
  default), `~ item` (GPCM), `~ step` (steps common to the items), `~ 1`
  (PCM); `~ item + step` (Rank-1, a<sub>il</sub> = α<sub>i</sub>
  γ<sub>l</sub>) is fitted with `xxirt_tppcm()`. Any other equality pattern
  goes through `index`.
* **Parameter table.** `irt_pars()` prints step discriminations and
  difficulties with standard errors and a model-fit header.
* **Model fit.** `m2()` gives limited-information fit statistics (M2, RMSEA,
  SRMSR, CFI, TLI) for any step design.
* **Step information.** `step_info()` splits item information into the
  contributions of each step.
* **Tests between models.** `wald_test()` tests equal steps within an item.
  `score_test()`, `item_test()` and `mi()` test a fitted model against a
  larger one without fitting it.
* **Rank-1.** `xxirt_tppcm(dat, design = ~ item + step)` fits
  a<sub>il</sub> = α<sub>i</sub> γ<sub>l</sub> and non-negativity bounds
  through `sirt::xxirt()`.
* **Differential step functioning.** On a multigroup fit,
  `dif_test(fit, g, by = "step")` screens every category boundary with a
  joint score (LM) test; `dsf()` decomposes flagged steps into CBD and
  location differences by likelihood ratio tests; `dsf_design()` builds
  configural, metric, scalar and partial multigroup models (lavaan-style
  `group.equal` / `group.partial`). `tppcmtree()` grows DIF trees when the
  grouping is unknown. Tables of tests report Holm-adjusted `p_adj`.
* **Components.** `get_parts()` extracts the quantities behind these
  functions (scores, information, Jacobian, posterior and more).
* **Simulation.** `sim_tppcm()` generates responses from a TPPCM.

## Workflow

Fit with TAM, then pass the fitted object to tppcm.

```r
library(tppcm)   # attaches TAM

set.seed(1)
a <- outer(c(0.8, 1.0, 1.2, 1.5, 1.0), c(1, 1.5, 0.7))
b <- matrix(c(-1, 0, 1), 5, 3, byrow = TRUE)
dat <- sim_tppcm(1000, disc = a, diff = b)

# 1. Fit the models with TAM
m_gpcm  <- tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
m_tppcm <- tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)

# 2. Parameters, fit and step-level tests
irt_pars(m_tppcm)
m2(m_tppcm)
wald_test(m_tppcm)

# 3. Where does the GPCM fail?
sp <- get_parts(m_gpcm)
score_test(sp)
item_test(sp)
mi(sp)

# 4. Rank-1: is there a step effect, an item x step interaction?
score_test(sp, against = ~ item + step)        # GPCM vs Rank-1
m_rank1 <- xxirt_tppcm(dat, design = ~ item + step)
irt_pars(m_rank1, restricted = TRUE)
score_test(m_rank1)                            # Rank-1 vs saturated

# 5. Differential step functioning between two groups
g <- rep(c("A", "B"), each = 500)
m_mg <- tam.mml.3pl(dat, E = tppcm(dat), group = g, est.variance = FALSE, verbose = FALSE)
st <- dif_test(m_mg, g, by = "step")           # screen every step (2 df)
flagged <- st$unit[st$p_adj < 0.05]
if (length(flagged)) dsf(m_mg, items = flagged) # CBD vs location of flagged steps
```

## Documentation

Two vignettes, each in English and in Traditional Chinese:

* `vignette("tppcm-tutorial", package = "tppcm")` (`"tppcm-tutorial-zh"`):
  the model and the package.
* `vignette("tppcm-recipes", package = "tppcm")` (`"tppcm-recipes-zh"`):
  short recipes for every function.

Cite with `citation("tppcm")`.

## References

* Yu, M.-N. (1991). *A two-parameter partial credit model* (Doctoral
  dissertation). University of Illinois at Urbana-Champaign.
* Masters, G. N. (1982). A Rasch model for partial credit scoring.
  *Psychometrika, 47*, 149-174.
* Muraki, E. (1992). A generalized partial credit model: Application of an
  EM algorithm. *Applied Psychological Measurement, 16*, 159-176.
