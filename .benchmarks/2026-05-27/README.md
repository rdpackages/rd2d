# rd2d Numerical and Speed Benchmarks

Generated workspace for comparing public `rd2d` releases against the local
working tree across R, Python, and Stata.

The intended invariant is:

- public `0.2.0`/`0.2.0.post1` is treated as legacy separate fitting;
- local `1.0.0` is called with `fitmethod = "separate"` / `fitmethod(separate)`
  for version-stability comparisons;
- local `fitmethod = "joint"` is benchmarked separately as the new default.

Generated CSVs and logs are intentionally kept under this workspace.
Generated directories (`inputs/`, `results/`, `profiles/`, `logs/`) and the
public Stata archive are ignored by the benchmark-local `.gitignore`. If this
repository is using a local rule that ignores `.gitignore` files, force-add
`.benchmarks/2026-05-27/.gitignore` when committing the benchmark scripts.

## How to rerun

From the repository root:

```powershell
python .benchmarks\2026-05-27\validate_workspace.py --require-results
python .benchmarks\2026-05-27\check_numerics.py --check-only
python .benchmarks\2026-05-27\check_numerics.py --skip-inputs
python .benchmarks\2026-05-27\check_numerics.py --skip-inputs --profiles
```

The first command validates the benchmark workspace and current generated
summaries. The second refreshes summaries and applies the tolerance guardrail
using existing benchmark outputs. The third regenerates the numerical benchmark
outputs. The fourth also runs the larger speed profiles.

The underlying individual commands are:

```powershell
Rscript .benchmarks\2026-05-27\make_inputs.R .
Rscript .benchmarks\2026-05-27\bench_r.R . public
Rscript .benchmarks\2026-05-27\bench_r.R . local_separate
Rscript .benchmarks\2026-05-27\bench_r.R . local_joint
python .benchmarks\2026-05-27\bench_python.py . public
python .benchmarks\2026-05-27\bench_python.py . local_separate
python .benchmarks\2026-05-27\bench_python.py . local_joint
Start-Process -FilePath 'C:\Program Files\Stata16\StataMP-64.exe' -ArgumentList @('/e','do','.benchmarks\2026-05-27\bench_stata.do',"`"$PWD`"",'public') -WorkingDirectory $PWD -WindowStyle Hidden -Wait
Start-Process -FilePath 'C:\Program Files\Stata16\StataMP-64.exe' -ArgumentList @('/e','do','.benchmarks\2026-05-27\bench_stata.do',"`"$PWD`"",'local_separate') -WorkingDirectory $PWD -WindowStyle Hidden -Wait
Start-Process -FilePath 'C:\Program Files\Stata16\StataMP-64.exe' -ArgumentList @('/e','do','.benchmarks\2026-05-27\bench_stata.do',"`"$PWD`"",'local_joint') -WorkingDirectory $PWD -WindowStyle Hidden -Wait
python .benchmarks\2026-05-27\summarize_benchmarks.py .
python .benchmarks\2026-05-27\check_tolerances.py .
python .benchmarks\2026-05-27\report_status.py
Rscript .benchmarks\2026-05-27\profile_r.R .
python .benchmarks\2026-05-27\profile_python.py . separate
python .benchmarks\2026-05-27\profile_python.py . joint
Start-Process -FilePath 'C:\Program Files\Stata16\StataMP-64.exe' -ArgumentList @('/e','do','.benchmarks\2026-05-27\profile_stata.do',"`"$PWD`"",'separate') -WorkingDirectory $PWD -WindowStyle Hidden -Wait
Start-Process -FilePath 'C:\Program Files\Stata16\StataMP-64.exe' -ArgumentList @('/e','do','.benchmarks\2026-05-27\profile_stata.do',"`"$PWD`"",'joint') -WorkingDirectory $PWD -WindowStyle Hidden -Wait
Start-Process -FilePath 'C:\Program Files\Stata16\StataMP-64.exe' -ArgumentList @('/e','do','.benchmarks\2026-05-27\test_stata_single_mlib.do',"`"$PWD\.benchmarks\2026-05-27\tmp_stata_single_mlib`"') -WorkingDirectory $PWD -WindowStyle Hidden -Wait
```

## Current baseline

The current run writes results under `results/` and Stata logs under `logs/`.

- Local `fitmethod = "separate"` is aligned across platforms, including the
  clustered/HC1 cases. Current max differences are:
  R vs Python `3.93e-07` absolute / `3.09e-11` relative; R vs Stata
  `2.08e-06` absolute / `1.64e-10` relative; Python vs Stata `1.11e-02`
  absolute / `2.64e-10` relative, with the absolute maximum coming from a very
  large covariance entry.
- Non-cluster public vs local `fitmethod = "separate"` remains exact for
  Python and Stata on common cells and stable for R up to `3.7e-12`.
- Public Python/Stata clustered cases now differ from local separate because
  local has been aligned to the R finite-sample clustered covariance formulas.
  Public cross-platform clustered differences are retained as legacy baseline
  evidence.
- Local `fitmethod = "joint"` differs from local `separate` where expected,
  mostly HC1/clustering and bandwidth-constant paths.

The benchmark scripts use Stata 16 by default. They pin public Stata to the
user-level `v0.2.0` ado install. The public Stata install needed the support
files from the `v0.2.0` tag (`rd2d_functions.do`, `lrd2d.mlib`, and package
companions) copied into `C:\Users\cattaneo\ado\plus\r`. Local development and
distribution now use a single `lrd2d.mlib` runtime artifact; the source
`rd2d_functions.do` is only for rebuilding that library.

## Speed Profile Notes

`profile_r.R`, `profile_python.py`, and `profile_stata.do` create larger
synthetic local workloads and write timings plus language-specific profile
summaries under `profiles/`. The first speed pass targeted clustered distance
covariance construction:

- R `dist_fuzzy_cluster_cov` median elapsed time dropped from about `1.28s`
  before optimization to the low tenths of a second after replacing repeated
  per-cluster row selection with shared cluster-sum aggregation, avoiding
  repeated treatment-side data-frame splits in location fits, and fast-pathing
  common kernel dispatch. R `loc_fuzzy_cluster_cov` is also in the low tenths
  of a second on the larger profile workload.
- R covariance-side `crossprod()` rewrites were retained after passing exact
  tests and the public/local benchmark guardrail. Coefficient-chain regrouping
  via `crossprod()` was tested but not retained because it moved one strict
  numerical-regression comparison beyond tolerance.
- A subsequent R data-management pass reduced repeated side splitting,
  bandwidth-row data-frame handling, outcome-matrix allocation, and temporary
  basis `cbind()` creation. Current larger-profile medians are about `0.09s`
  for `loc_sharp_cov`, `0.17s` for `loc_fuzzy_cluster_cov`, `0.16s` for
  `loc_bw`, `0.05s` for `dist_sharp_cov`, `0.14s` for
  `dist_fuzzy_cluster_cov`, and `0.09s` for `dist_bw`.
- Python and Stata now use the same cluster-sum strategy in their local
  clustered covariance paths, and both reuse full cluster groups across
  repeated evaluation fits.
- Python additionally reuses full cluster groups inside repeated evaluation
  fits and uses a NumPy fast path for cluster coding. On the larger profile
  workload, Python median times are about `0.043s`/`0.052s` for separate
  location/distance fuzzy clustered covariance and `0.042s`/`0.056s` for the
  joint counterparts.
- The Stata larger workload is now substantially faster in clustered
  covariance cases than the original baseline: separate location/distance
  fuzzy clustered covariance medians are about `0.168s`/`0.212s`; joint
  medians are about `0.170s`/`0.220s` under Stata 16. A side-specific
  clustered distance-fit branch was retained after passing the numerical
  guardrail, and the same side-specific branch is now used for small
  noncluster distance batches. The small noncluster location path now also
  uses side-specific fits when covariance adjustment and clustering are absent,
  which removes the largest synthetic location fixed-bandwidth gap while
  preserving public/local and cross-platform numerical agreement. Larger
  noncluster distance batches keep the shared-design path because the
  side-specific branch did not improve the larger profile workload there.
  Location bandwidth selection now caches repeated side splits, side-specific
  cluster vectors, weight-count calculations, and non-covariate side outcome
  vectors; current separate/joint `loc_bw` medians are about `0.199s` under
  Stata 16. Local fits skip leverage calculations unless HC2/HC3 actually
  needs them. Distance fits use a nonnegative-kernel fast path for
  side-specific distance weights.
- After the speed changes, the local `fitmethod = "separate"` numerical
  comparison remains stable at the same cross-platform tolerances reported
  above.

## Tolerance Guardrail

`check_tolerances.py` refreshes the summary CSVs and fails if either of the
core numerical invariants is violated:

- local `fitmethod = "separate"` must agree across R, Python, and Stata within
  `1e-8` relative tolerance;
- non-cluster public-vs-local `fitmethod = "separate"` comparisons must agree
  within `1e-8` absolute and relative tolerance.

Clustered public-vs-local Python/Stata comparisons are intentionally excluded
because local `1.0.0` now uses the R-aligned finite-sample clustered covariance
formulas.

## Compact Status Report

`report_status.py` writes `results/status_report.md` and
`results/status_report.json` from the current summary and profile CSVs. These
reports are generated artifacts, but they are useful for a quick audit after
running `check_numerics.py`.
