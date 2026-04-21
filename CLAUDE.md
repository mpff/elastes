# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package overview

`elastes` is an R package for computing **elastic full Procrustes means** of open planar curves that are only sparsely and irregularly observed. The mathematical approach treats planar curves as complex-valued functions and estimates a Procrustes mean from a smoothed Hermitian covariance surface. See Stöcker et al. (2022, arXiv:2203.10522) and Steyer, Stöcker & Greven (2022, doi:10.1111/biom.13706).

The package is standard R-package layout (DESCRIPTION/NAMESPACE/R/man/tests/vignettes) with `roxygen2` used to generate `man/*.Rd` and `NAMESPACE` from roxygen comments in `R/*.R`. Do not edit `NAMESPACE` or `man/*.Rd` by hand.

## Common commands

All commands run from the package root. These rely on `devtools`/`testthat`/`roxygen2`/`rcmdcheck`:

```r
# Regenerate NAMESPACE and man/ pages after editing roxygen comments
devtools::document()

# Full test suite (testthat edition 3)
devtools::test()

# Run a single test file
devtools::test(filter = "compute_elastic_shape_mean")   # matches test-compute_elastic_shape_mean.R

# Load package for interactive work (like library(elastes) but from source)
devtools::load_all()

# Rebuild README.md from README.Rmd (README.md is generated — never edit directly)
devtools::build_readme()

# CRAN-style check locally (matches GitHub Actions R-CMD-check workflow)
rcmdcheck::rcmdcheck(args = c("--no-manual", "--as-cran"))

# Coverage (used by codecov workflow)
covr::package_coverage()
```

CI runs R-CMD-check on Windows + Ubuntu (devel/release/oldrel-1) via `.github/workflows/R-CMD-check.yaml`. `pkgdown` site and `test-coverage` are separate workflows.

## Architecture

### Exported API (see NAMESPACE)

Only three user-facing functions are exported — everything else in `R/` is internal:

- `compute_elastic_shape_mean()` — main entry point; constructor for class `elastic_shape_mean`.
- `get_procrustes_fit()` — recover the Procrustes-aligned version of an original data curve from the mean object.
- `get_evals()` — S3 generic with methods for `default`, `data.frame`, and `elastic_shape_mean` (evaluates a curve or the mean on a grid).
- `plot.elastic_shape_mean()` — S3 plot method.

### Estimation pipeline

`compute_elastic_shape_mean()` (`R/compute_elastic_shape_mean.R`) is orchestration: input validation, arc-length parameterisation, centering, unit-length prescaling, conversion to SRV (square-root-velocity) representation via `elasdics::get_srv_from_points`. It then delegates the actual iterative estimation to `fit_mean()`.

`fit_mean()` (`R/fit_mean.R`) is the core algorithm. Each iteration:

1. `get_model_data_complex()` — build complex-valued SRV observations on the current warped parametrisation (differs for `type = "smooth"` vs `"polygon"`).
2. Build pairwise product data `qq = Conj(q_i) * q_j` on an (s,t) grid per curve.
3. `smooth_cov_surface()` — tensor-product P-spline smoothing of the complex covariance surface via `mgcv::bam`. Real part uses a symmetric `bs="symm"` basis; imaginary part uses a skew-symmetric one. Measurement-error variance is controlled by `var_type` (`"smooth"`, `"constant"`, `"zero"`).
4. `get_complex_coef_matrix()` → eigendecomposition of `L.inv' · beta · L.inv` (where `L = chol(G^-1)`, `G` = Gram matrix of mean basis) gives the Procrustes mean coefficients as the leading eigenvector.
5. `fit_alignment_proc2d()` — per-curve optimal rotation/scaling alignment. The `"smooth"` branch uses a conditional-score formulation in the covariance eigen-basis; the `"polygon"` branch does an interval integration against the current mean function.
6. `get_optimal_t()` — warping step, dispatches to `elasdics::find_optimal_t` (smooth) or `find_optimal_t_discrete` (polygon).
7. Convergence when relative L2 change of mean coefficients `< eps`.

### Basis machinery

`R/basis.R` provides `get_knots()`, `make_design()` (via `splines::splineDesign`), and `get_gram_matrix()` (uses `orthogonalsplinebasis` for smooth/order-2 B-splines, a diagonal `diag(diff(knots))` for the polygon/piecewise-constant case). `type` is the single switch between smooth (linear SRV splines ⇒ differentiable mean) and polygon (piecewise constant SRV ⇒ piecewise linear mean); spline degree is `ifelse(type == "smooth", 1, 0)`.

### Dependency on `elasdics`

Several internal functions from the `elasdics` package are re-bound in `R/elasdics_imports.R` via `utils::getFromNamespace`: `remove_duplicate`, `srvf_to_curve`, `find_optimal_t`, `find_optimal_t_discrete`. Public `elasdics` functions (`get_arc_length_param`, `get_srv_from_points`, `get_points_from_srv`, `center_curve`) are called with `elasdics::` prefix. If bumping the `elasdics` version, check both call sites.

### Result object

`compute_elastic_shape_mean()` returns an S3 `elastic_shape_mean` list with `coefs`, `knots`, `type`, `variance`, per-curve `data_curves` (each carrying `rotation`, `scale`, `translation`, `polygon_length`, `norm_factor`, `dist_to_mean` as attributes), and a `fit` sublist with diagnostics (`cov_fit`, `cov_pca`, `pfit_coefs`, `G_optims`, `b_optims`, `l_optims`, `n_iter`, `gram`). Attributes on `data_curves[[j]]` are what `get_procrustes_fit()` reads to invert translation/rotation/scale.

## Conventions

- Roxygen `@export` tags drive `NAMESPACE` — add or remove them rather than editing `NAMESPACE` directly, then run `devtools::document()`.
- `README.md` is generated from `README.Rmd`; regenerate with `devtools::build_readme()`, do not edit `README.md` by hand.
- Tests live under `tests/testthat/` and follow `test-<topic>.R` naming; the package uses `testthat` edition 3 (per `Config/testthat/edition` in DESCRIPTION).
- Several parameters (`pfit_method`, `var_type`, `smooth_warp`, `cluster`) are marked `(experimental)` in the docs — preserve that label when touching their documentation.
