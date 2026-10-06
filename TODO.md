# tppcm: open items (as of 2026-10-06)

## Interface / consistency
- [ ] Unify `irt_pars(IRTpars = )` with `param = c("irt", "si")` used elsewhere.
- [ ] `score_test()`, `item_test()`, `mi()`: add `info = "robust"` (currently louis/opg only; `dif_test()` defaults to robust). Consider merging them into `score_test(by = c("model", "item", "step"))`; `mi()` gives duplicate rows for two-step items.
- [ ] `dif_test()`: message when a numeric covariate is the grouping of a multigroup fit (silently gives CvM instead of the split LM).
- [ ] `dif_test()`: trim `parm` aliases (`"all"`, `"slope"`/`"int"` under `param = "irt"`)?
- [ ] `xxirt_vcov()`: turn into a `vcov.xxirt_tppcm` method (sirt vcov workaround).
- [ ] `dsf_design(group.equal = "ibd")` (IBDs equal, CBDs free) is accepted but undocumented, and `dsf()` refuses it.
- [ ] `m2()`: M2 is sensitive to the convergence criterion (changes of 12-40 on the SEL data with stricter conv while logLik changes < 0.07); warn or document.
- [ ] One test in test-dif.R prints the impact message to the console.

## DSF
- [ ] Purification / anchors: `anchor =` (free baseline) for `dsf()`; iterative purification (Khalid & Glas, 2014).
- [ ] DSF effect sizes (area between step curves, expected-score differences).
- [ ] Simulate more than two groups and misfit present in one group only.
- [ ] Optional MH-type DSF (Penfield, 2007) matched on the TPPCM weighted score (prototype: paper_B_zh/draft_v4/code/mfr/mh1.R).

## Model selection (G-DINA analogy)
- [ ] `select_design()`: item-wise choice of the CBD structure from a saturated fit (flat vs free) by Wald (or LR) tests, in the spirit of de la Torre & Lee (2013) / GDINA::modelcomp(); returns `index` and the refitted mixed model.
- [ ] Item-wise "common shape" (Rank-1) option: nonlinear Wald, mixed models via xxirt.
- [ ] Note the G-DINA analogy in docs: the saturated model is link-invariant, main-effects-only reductions depend on the link (A-CDM / LLM / R-RUM vs additive / Rank-1).
- [ ] Explanatory step discriminations (LLTM-type): item covariates in `design` (e.g. `~ step * reversed`).
- [ ] Verify references before use: de la Torre (2011), de la Torre & Lee (2013), Ma, Iaconangelo & de la Torre (2016), Sorrel et al. (2017).

## Infrastructure
- [ ] PDF manual not built locally (TeX lacks inconsolata.sty).
