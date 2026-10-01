# tppcm

Step discriminations for the two-parameter partial credit model (TPPCM),
built on top of **TAM**.

The TPPCM was proposed by Yu (1991): every step *l* of item *i* in Masters'
(1982) partial credit model gets its own discrimination *a<sub>il</sub>*:

P(X<sub>i</sub> = k | θ) ∝ exp( Σ<sub>l ≤ k</sub> (a<sub>il</sub> θ − d<sub>il</sub>) ),  k = 0, …, K−1.

The GPCM (Muraki, 1992) is the special case with equal steps within items.

The saturated TPPCM is a reparameterization of the nominal response model, so
this package does **not** add another estimator. Estimation stays in TAM
(and in sirt where TAM cannot express a model exactly). The package adds what
these packages do not provide.

| Task | Existing software | tppcm |
|---|---|---|
| Fit the PCM, GPCM | `tam.mml()`, `tam.mml.2pl()` | — |
| Fit the TPPCM with step discriminations as parameters; equal or fixed steps | `tam.mml.3pl()` needs a design array | `tppcm()` (alias `step_design()`) |
| Standard errors of step discriminations | `se.gammaslope` ignores the covariance with thresholds and is too small | `irt_pars()` (`long = TRUE` for intervals) |
| Compact item/group parameter table (mirt-style) with fit header | `summary()` is long | `irt_pars()` |
| Step contributions to item information (item and test information are in TAM) | `IRT.informationCurves()` gives item/test information only | `step_info()` |
| Limited-information fit: M2, RMSEA2, SRMSR, CFI, TLI | AIC/BIC, SRMSR in TAM; M2 only in mirt, not for step designs | `m2()` |
| Wald test: equal steps within an item | — | `wald_test()` |
| Score tests, modification indices, EPCs (no need to fit the larger model) | — | `score_test()`, `mi()`, `item_test()` |
| Product form a<sub>il</sub> = α<sub>i</sub> γ<sub>l</sub>; bound a<sub>il</sub> ≥ 0 | not expressible in TAM | `xxirt_tppcm()` (via `sirt::xxirt()`) |
| Score-based DIF tests per step, item or block, discriminations and difficulties separately | — for TAM (psychotools / mirt: other models) | `dif_test()` |
| DIF trees on step discriminations or on locations separately | `psychotree::gpcmtree()` (item-level slopes, slope-intercept form) | `tppcmtree(parm = "disc" / "diff")` |
| Covariance of `xxirt` fits with shared parameters | `vcov()` fails in sirt 4.2-133 | `xxirt_vcov()` |

## Documentation

Two vignettes, each in English and in Traditional Chinese:

* `vignette("tppcm-tutorial", package = "tppcm")` (`"tppcm-tutorial-zh"`):
  the model, the package and the reasoning behind each tool.
* `vignette("tppcm-recipes", package = "tppcm")` (`"tppcm-recipes-zh"`):
  recipes (problem, solution, discussion, see also) that use every function,
  including deliberately wrong inputs and the messages they produce.

The examples use data sets of other packages: `CDM::data.Students`
(mathematics self-concept, 0–3, with covariates), `TAM::data.timssAusTwn.scored`
(TIMSS 2011, mixed dichotomous and partial-credit items, two countries),
`psychotools::VerbalAggression` (24 items, 0–2, gender and anger),
`mirt::Science` and `TAM::data.gpcm`.

Cite with `citation("tppcm")`.

## Using TAM's own tools on the fit

In TAM 4.3-25, `tam.wle()`, `tam.fit()` and `plot()` fail on any
`tam.mml.3pl` object (with or without a design array) with
`object 'res' not found`, because they dispatch only for class `tam.mml`.
`tam.modelfit()`, `IRT.compareModels()` and `IRT.informationCurves()` work.
Two workarounds:

```r
tam.mml.wle2(fit)                                  # WLE person estimates
class(fit) <- c("tam.mml.3pl", "tam.mml"); tam.wle(fit); tam.fit(fit); plot(fit)
```

## Installation

```r
devtools::install("tppcm", build_vignettes = TRUE)   # from the directory containing the package
```

## Parameter names

Outputs use words, not letters, because `d` means different things in
different packages:

| Form | tppcm | Elsewhere |
|---|---|---|
| `disc * (theta - diff)` | `disc` (step discrimination), `diff` (step difficulty) | Masters' step difficulty; TAM `beta + tau` (PCM `xsi`); mirt `b` (`IRTpars = TRUE`); psychotools `threshpar()` |
| `slope * theta + int` | `slope`, `int` (step intercept) | cumulated over steps: mirt / psychotools category intercepts `d_k` (= `-AXsi_` in TAM) |

## Models

`tppcm(dat, model = )` gives the common designs; `index` gives any linear
design (same index = equal step discriminations).

```
            tppcm               a_il free  (default)
             |
            rank1               a_il = alpha_i gamma_l   (xxirt_tppcm(dat, "rank1"))
          /    \
  gpcm (alpha_i)  step (gamma_l)
          \    /
            pcm                 a_il = a  (= PCM with free variance)
```

## Design: fitted object in, components out

The package never replaces the TAM (or xxirt) object. Every function takes
the fitted object directly; the responses are read from the fit.

```r
fit <- tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE)
irt_pars(fit); wald_test(fit); m2(fit)      # works on the TAM object
```

`get_parts(fit, what)` extracts the mathematical components behind these
functions, like `lme4::getME()`, `lavaan::lavInspect()` or
`mirt::extract.mirt()`. Without `what` it returns all of them; pass that
object instead of `fit` to avoid recomputing. It has `vcov()`, `logLik()`
and `nobs()` methods.

```r
get_parts(fit, "disc")              # step discriminations
get_parts(fit, c("loglik", "npar"))
sp <- get_parts(fit)                # all components, computed once
score_test(sp); mi(sp)
```

| `what` | Content |
|---|---|
| `"disc"`, `"diff"`, `"int"` | items × steps matrices of discriminations, difficulties, intercepts |
| `"est"`, `"se"`, `"vcov"` | slope parameters of the fitted model (item slopes, α, γ, …), SEs, covariance |
| `"par"`, `"partable"` | saturated parameters (per item a<sub>1..L</sub>, d<sub>1..L</sub>) and their description |
| `"vcov.sat"`, `"se.sat"` | implied covariance and SEs of the saturated parameters |
| `"information"`, `"information.opg"` | observed (Louis) and outer-product information of the saturated model |
| `"gradient"`, `"scores"` | score vector and case-wise scores of the saturated model |
| `"jacobian"`, `"curvature"` | Jacobian of the fitted model in the saturated space, curvature term |
| `"posterior"`, `"nodes"`, `"prior"`, `"probs"` | posterior, quadrature nodes, latent distribution, response probabilities |
| `"data"`, `"nobs"`, `"loglik"`, `"npar"` | responses, persons, log-likelihood, parameters |
| `"model"`, `"groups"`, `"design"` | model type, latent means and variances, design array `E` |

The model is detected from the fit:

| Fit | Model |
|---|---|
| `tam.mml.3pl(E = tppcm(...))` | linear step design (`"E"`), fixed steps detected |
| `tam.mml.2pl(irtmodel = "2PL")` | saturated TPPCM (`"tppcm"`) |
| `tam.mml.2pl(irtmodel = "GPCM")` | GPCM (`"gpcm"`) |
| `tam.mml()` | PCM with free latent variance (`"pcm"`) |
| `xxirt_tppcm()` | `"tppcm"`, `"gpcm"`, product form `"rank1"` |

Also accepted: binary items, `tam.mml(irtmodel = "RSM")` and other location
designs, several groups, person weights (`pweights`), latent regression (with
a warning). Refused with a message: guessing parameters, `tam.mml.mfr()`,
`tam.jml()`, multidimensional fits, objects from other packages.

## Workflow

```r
library(tppcm)   # attaches TAM

set.seed(1)
a <- outer(c(0.8, 1.0, 1.2, 1.5, 1.0), c(1, 1.5, 0.7))    # true product form
b <- matrix(c(-1, 0, 1), 5, 3, byrow = TRUE)
dat <- sim_tppcm(1000, disc = a, diff = b)

# 1. Fit with TAM
m_pcm   <- tam.mml.3pl(dat, E = tppcm(dat, "pcm"), est.variance = FALSE, verbose = FALSE)  # PCM, common slope
m_gpcm  <- tam.mml.2pl(dat, irtmodel = "GPCM", verbose = FALSE)
m_tppcm <- tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE, verbose = FALSE)
m_step  <- tam.mml.3pl(dat, E = tppcm(dat, "step"), est.variance = FALSE, verbose = FALSE)  # a_il = gamma_l
IRT.compareModels(m_pcm, m_gpcm, m_tppcm)

# 2. Clean parameter table and step discriminations with correct standard errors
irt_pars(m_tppcm)                 # fit (AIC, BIC, M2, RMSEA, SRMSR, CFI, TLI), disc, diff, SEs, groups
m2(m_tppcm)                       # limited-information fit alone
irt_pars(m_tppcm, long = TRUE)    # long format with confidence intervals
wald_test(m_tppcm)                     # per item: a_i1 = a_i2 = a_i3 ?

# 3. Start from the GPCM: where does it fail?
sp <- get_parts(m_gpcm)           # components, computed once
score_test(sp)                  # GPCM vs saturated TPPCM
score_test(sp, against = "rank1")  # GPCM vs product form
item_test(sp)                  # which items
mi(sp)                          # which steps (MI and EPC)

# 4. Product form (sirt::xxirt; same nodes as TAM, so log-likelihoods compare)
m_rank1 <- xxirt_tppcm(dat, "rank1")
irt_pars(m_rank1, restricted = TRUE)    # alpha, gamma with SEs
score_test(m_rank1)                  # product form vs saturated
IRT.compareModels(m_gpcm, m_rank1, m_tppcm)$LRtest   # use LRtest; not anova() with mixed objects

# Restrictions on steps in TAM: same index = equal; gammaslope.fixed = fixed
idx <- rbind(1:3, 4:6, rep(7, 3), rep(8, 3), rep(9, 3))   # items 3-5 GPCM
m_mix <- tam.mml.3pl(dat, E = tppcm(dat, index = idx), est.variance = FALSE, verbose = FALSE)
wald_test(m_mix)                       # only items 1-2 are tested
```

## Model fit

`irt_pars()` starts with a fit header: log-likelihood, number of parameters,
AIC and BIC (as `IRT.IC()` in TAM), and the output of `m2()`: M2 with its
p-value, RMSEA2 with a 90% interval, SRMSR, CFI and TLI (independence model as
baseline). `m2()` differentiates the margins with respect to the saturated
parameters and maps them with the Jacobian of the fit, so it works for every
step design and the product form. Types: `"M2"` (all univariate and bivariate
margins; default up to 2000 margins), `"C2"`, `"M2*"`.

Checks: identical to `mirt::M2()` for the dichotomous 2PL; `"C2"` correlates
.9996 with mirt for the GPCM. Under a true GPCM (8 items, 4 categories,
N = 1000, 100 replications) the mean statistic and rejection rate were
243.98 / .05 for M2 (df 244), 19.90 / .04 for C2 (df 20) and 4.21 / .05 for
M2* (df 4). Missing responses are allowed: margins use the available cases
and their covariance is scaled by the joint numbers of observed persons
(valid under MCAR; with 15% missing the means were 244.98, 19.86 and 3.88).

The reduced statistics (`"C2"`, `"M2*"`) collapse categories into item-score
moments and can miss unequal step discriminations; use the full M2 for TPPCM
questions.

Note: in TAM fits `pi.k` (the `prob.theta` attribute of `IRT.posterior()`)
is the mean posterior, not the normal prior of the likelihood; `m2()` uses the
fitted normal distribution.

## Score-based DIF tests

`dif_test(fit, covariate)` tests whether parameters of a single-group fit
vary with a person covariate (Merkle & Zeileis, 2013), without fitting a
multiple-group model. Units are single parameters (`by = "param"`), items or
the whole block; each unit uses efficient scores (residualized on all other
parameters), so DIF in the difficulties does not leak into tests of the
discriminations and vice versa. Factors give the LM test, ordered factors
`"maxLMo"`, numeric covariates `"CvM"`; p-values are also Holm-adjusted.

```r
fit  <- tam.mml.3pl(dat, E = tppcm(dat), est.variance = FALSE)
sp <- get_parts(fit)
dif_test(sp, group, parm = "disc")               # every step discrimination
dif_test(sp, group, parm = "diff", by = "item")  # difficulties per item
```

**Impact is not DIF.** A single-group fit assumes one latent distribution.
When the groups differ in ability (or a numeric covariate is correlated with
it), the difficulty scores of all items move with the covariate: with a 0.5
SD mean difference and no DIF, the difficulty tests rejected in 98–100% of
50 simulated samples (the discrimination tests in 4%). Fit the model with
`group = covariate` (or `Y = covariate`) and run `dif_test()` on that fit;
the group means and variances are then part of the model and the false
alarms vanish (0% difficulty, 12% discrimination with mean and variance
impact) while the power for real step DIF is retained. `dif_test()` prints a
reminder when a factor is tested on a single-group fit.

Simulation (6 items, 4 categories, N = 1500, binary covariate, 200
replications; rate of at least one Holm-adjusted rejection):

| DIF in the covariate | TPPCM fit: steps, disc | TPPCM: items, disc | TPPCM: items, diff | GPCM fit: items, disc | GPCM: items, diff |
|---|---|---|---|---|---|
| none | .07 | .02 | .05 | .04 | .05 |
| one step flatter | 1.00 | 1.00 | .05 | .98 | .41 |
| two steps in opposite directions (item mean unchanged) | .77 | .84 | .06 | .07 | .20 |
| whole item discrimination | .87 | 1.00 | .06 | 1.00 | .05 |
| difficulty | .04 | .03 | 1.00 | .03 | 1.00 |

With a GPCM fit, step-discrimination DIF is either invisible or reported as
difficulty DIF; fit the saturated TPPCM to test it.

## DIF trees

`tppcmtree()` grows a model-based recursive partitioning tree
(`partykit::mob()`) with a TAM step-design fit in every node. `parm = "disc"`
splits only on instability of the discriminations, `parm = "diff"` only on
instability of the step difficulties (`"slope"` / `"int"` with
`param = "si"`). Node models: `"tppcm"`, `"gpcm"`, `"pcm"`. Testing both together lets a strong
location DIF mask a step-discrimination DIF.

```r
df$resp <- resp_matrix
tppcmtree(resp ~ gender + age, data = df, parm = "disc")  # discrimination DIF
tppcmtree(resp ~ gender + age, data = df, parm = "diff")  # difficulty DIF
tppcmtree(resp ~ gender + age, data = df, parm = "diff", impact = "gender")
get_parts(fit, "estfun")                                # the case-wise scores used
```

Trees based on marginal maximum likelihood split on a covariate related to
ability even without DIF (with a 0.5 SD mean impact, `pcm`/`diff` and
`gpcm`/`all` trees split in 100% of simulated samples; CML trees such as
`pctree` do not). `impact = ` names a grouping factor whose latent means and
variances are estimated in every node (`group = ` in TAM), as in
`psychotools::gpcmodel(impact = )`, so that splits on it reflect item
differences only.

By default (`orthogonal = TRUE`) the tested block of scores is replaced by
efficient scores (residualized on the other block); otherwise `mob()`'s
jointly decorrelated process lets location DIF leak into the slope test.
Simulation (5 items, N = 2000, 60 replications, root-node rejection rates):

| | discrimination-DIF covariate | difficulty-DIF covariate | noise |
|---|---|---|---|
| `parm = "disc"`, not orthogonal | .97 | .55 | .02 |
| `parm = "disc"`, orthogonal | .95 | .00 | .02 |
| `parm = "diff"`, not orthogonal | .07 | 1.00 | .02 |
| `parm = "diff"`, orthogonal | .03 | 1.00 | .02 |

Without DIF (6 items, N = 1200, 100 replications) the splitting rates were
.04 (`"disc"`), .05 (`"diff"`) and .09 (`"all"`).

Comparison with psychotree (6 items, 4 categories, N = 1500, binary covariate
`z` plus a noise covariate, 40 replications; rate of a root split on `z`; the
noise split rate was at most .075 for all methods):

| Data / DIF in `z` | `pctree` | `tppcmtree` pcm, `diff` | `gpcmtree` | `tppcmtree` gpcm, `all` | `gpcmtree` slopes only | `tppcmtree` gpcm, `disc` | `tppcmtree` gpcm, `diff` | `tppcmtree` tppcm, `disc` | `tppcmtree` tppcm, `diff` |
|---|---|---|---|---|---|---|---|---|---|
| PCM, none | .08 | .08 | .05 | .05 | .03 | .05 | .08 | .00 | .03 |
| PCM, difficulty | 1 | 1 | 1 | 1 | **.25** | **.00** | 1 | .00 | 1 |
| GPCM, none | .10 | .10 | .08 | .08 | .03 | .00 | .10 | .00 | .05 |
| GPCM, difficulty | 1 | 1 | 1 | 1 | **.20** | **.05** | 1 | .00 | 1 |
| GPCM, item discrimination | .23 | .05 | .93 | .93 | .83 | .97 | .03 | .85 | .05 |
| TPPCM, none | .03 | .00 | .00 | .00 | .00 | .05 | .05 | .00 | .03 |
| TPPCM, one step discrimination | .08 | .03 | .62 | .62 | .93 | .90 | .13 | .88 | .00 |
| TPPCM, two steps in opposite directions (item mean unchanged) | .10 | .10 | .10 | .10 | .05 | .03 | .10 | **.55** | .05 |
| GPCM, two steps in opposite directions (item mean unchanged) | .20 | .30 | .42 | .42 | .35 | .15 | .28 | **.75** | .05 |

* `tppcmtree(model = "pcm", parm = "diff")` matches `pctree` (CML) when the
  covariates are unrelated to ability (see the impact note above), and
  `tppcmtree(model = "gpcm", parm = "all")` matches `gpcmtree` exactly.
* Testing only the slopes of `gpcmtree` (`parm =` slope indices) lets
  difficulty DIF trigger slope splits (.20–.25); the orthogonalized
  `parm = "disc"` does not (.00–.05).
* Testing slopes alone gives more power for discrimination DIF than testing
  all parameters (.90 vs .62).
* `pctree` partly reads item-discrimination DIF as difficulty DIF (.23);
  the PCM node with `parm = "diff"` does not (.05).
* When two steps of an item change in opposite directions, the TPPCM node
  finds it best (.55, .75). PCM and GPCM trees either miss it (.03–.10) or
  split partly through the difficulties (.20–.42), i.e. they report it as
  location DIF.

## Methods in brief

* **Score core.** Person scores of the saturated parameters are posterior
  expectations of complete-data scores (Fisher's identity); the observed
  information is obtained with Louis' formula. Only `IRT.posterior()` and
  `IRT.irfprob()` of the fitted model are needed.
* **Restricted models** are maps *h(ϑ)* into the saturated space with
  Jacobian *J*. Their observed information is *J'IJ − Σ g<sub>k</sub> ∇²h<sub>k</sub>*;
  the second term matters only for the bilinear product form.
* **Score test** for extra directions *D*:
  *S = g<sub>D</sub>' (D'ID − D'IJ (J'IJ)<sup>−1</sup> J'ID)<sup>−1</sup> g<sub>D</sub>*.
  Against the saturated model this is *g'I<sup>−1</sup>g*. A modification index
  is the 1-df case for a single *a<sub>il</sub>*.

## Limitations

* The latent trait distribution is treated as fixed (N(0, 1) for the
  TPPCM, GPCM and product form; the PCM's variance is handled as a common
  slope).
* Items must have the same number of categories for `xxirt_tppcm()` and the
  product form; `tppcm()` allows different numbers.
* For likelihood-ratio comparisons, fit the PCM as `tppcm(dat, "pcm")` in
  `tam.mml.3pl()` (common slope, variance 1) rather than with `tam.mml()`
  (free variance). The two are the same model, but on TAM's default 21
  nodes their quadrature errors differ; with 16 items the log-likelihoods
  differed by up to 5 units and LR statistics against slope-based models
  became negative. Alternatively use more nodes.
* With estimates on the bound *a<sub>il</sub> = 0* the χ² reference
  distribution of LR and score tests is not valid.
* In `tam.mml.3pl()`, do not use `gammaslope.constr.V` for equality
  constraints (does not reach the MLE) or `userfct.gammaslope` for bounds
  (truncation is not the constrained MLE).
* In small samples (about fewer than 60 persons per item category), tests
  against the saturated model can be liberal; tests between adjacent models
  (e.g. GPCM vs product form) hold their level. (With 5 items, 4 categories
  and N = 250 the GPCM-vs-TPPCM score test rejected in 2.5% of 40 samples,
  so the effect is small.)
* `m2()` with missing responses assumes MCAR. When non-response depends on
  ability or on the responses (omits, not-reached items) the statistic is
  invalid: with missingness depending on one item's response, the true model
  was rejected in every simulated sample.
* The EPC of `mi()` is a one-step approximation and can overshoot the actual
  change when the MI is large (0.95 vs 0.64 in one check).
* Under misfit `m2(type = "C2")` can differ from mirt's C2 by 10–20% because
  the covariance matrix of the margins is taken from the model here and
  partly from the sample in mirt; under the true model the two agree.

## References

* Yu, M.-N. (1991). *A two-parameter partial credit model* (Doctoral
  dissertation). University of Illinois at Urbana-Champaign.
* Masters, G. N. (1982). A Rasch model for partial credit scoring.
  *Psychometrika, 47*, 149-174.
* Muraki, E. (1992). A generalized partial credit model: Application of an
  EM algorithm. *Applied Psychological Measurement, 16*, 159-176.
* Bock, R. D. (1972). Estimating item parameters and latent ability when
  responses are scored in two or more nominal categories. *Psychometrika,
  37*, 29-51.
* Maydeu-Olivares, A., & Joe, H. (2006). Limited information goodness-of-fit
  testing in multidimensional contingency tables. *Psychometrika, 71*,
  713-732.
* Merkle, E. C., & Zeileis, A. (2013). Tests of measurement invariance
  without subgroups. *Psychometrika, 78*, 59-82.
* Zeileis, A., Hothorn, T., & Hornik, K. (2008). Model-based recursive
  partitioning. *Journal of Computational and Graphical Statistics, 17*,
  492-514.
