# tppcm

<!-- badges: start -->
[![R-CMD-check](https://github.com/jiewenTsai/tppcm/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jiewenTsai/tppcm/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Step discriminations for the two-parameter partial credit model (TPPCM;
Yu, 1991), built on top of **TAM**.

In the TPPCM every step *l* of item *i* has its own discrimination
*a<sub>il</sub>*:

P(X<sub>i</sub> = k | θ) ∝ exp( Σ<sub>l ≤ k</sub> (a<sub>il</sub> θ − d<sub>il</sub>) ),  k = 0, …, K−1.

The generalized partial credit model (GPCM; Muraki, 1992) is the special case
with equal steps within items, and the partial credit model (PCM; Masters,
1982) has one common discrimination.

Estimation stays in TAM (and in sirt for the product form). tppcm takes the
fitted object and adds standard errors, tests, fit statistics and DIF tools
for step discriminations.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jiewenTsai/tppcm", build_vignettes = TRUE)
```

## Features

* **Step designs for TAM.** `tppcm()` builds the design array for
  `tam.mml.3pl()`: saturated TPPCM, GPCM, steps-only, PCM, or any
  equality pattern via `index`.
* **Parameter table.** `irt_pars()` prints step discriminations and
  difficulties with standard errors and a model-fit header.
* **Model fit.** `m2()` gives limited-information fit statistics (M2, RMSEA,
  SRMSR, CFI, TLI) for any step design.
* **Step information.** `step_info()` splits item information into the
  contributions of each step.
* **Tests between models.** `wald_test()` tests equal steps within an item.
  `score_test()`, `item_test()` and `mi()` test a fitted model against a
  larger one without fitting it.
* **Product form.** `xxirt_tppcm()` fits a<sub>il</sub> = α<sub>i</sub>
  γ<sub>l</sub> and non-negativity bounds through `sirt::xxirt()`.
* **DIF.** `dif_test()` gives score-based DIF tests per step, item or block,
  and `tppcmtree()` grows model-based DIF trees. Both treat discriminations
  and difficulties separately.
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

# 4. Product form
m_rank1 <- xxirt_tppcm(dat, "rank1")
irt_pars(m_rank1, restricted = TRUE)

# 5. DIF on step discriminations
group <- rep(c("A", "B"), each = 500)
dif_test(get_parts(m_tppcm), group, parm = "disc")
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
