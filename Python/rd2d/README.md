# rd2d

`rd2d` provides local polynomial estimation, robust bias-corrected inference,
and bandwidth helpers for boundary discontinuity designs with bivariate running
variables. The package includes location-based and distance-based methods,
sharp and fuzzy designs, pointwise confidence intervals, covariance-backed
summary inference, uniform confidence bands, and aggregate boundary effect
summaries.

```python
from rd2d import rd2d, rdbw2d, rd2d_dist, rdbw2d_dist, summary
```

## Main Functions

- `rd2d()`: location-based estimation and inference.
- `rdbw2d()`: location-based bandwidth selection.
- `rd2d_dist()` and `rd2d_distance()`: distance-based estimation and inference.
- `rdbw2d_dist()` and `rdbw2d_distance()`: distance-based bandwidth selection.
- `summary()`: summary tables with optional uniform bands, WBATE, and LBATE.

In the source repository, the sibling scripts `../rd2d_illustration.py` and
`../rd2d_plot.py` read `../rd2d_data.csv` and run the same dataset-based
illustration workflow as the R and Stata examples. They do not create an
`output/` folder or save fitted objects, tables, or plot files.

## Installation

```sh
python -m pip install rd2d
```

For local development:

```sh
python -m pip install -e .
```

Optional plotting and testing dependencies can be installed with:

```sh
python -m pip install -e ".[plots,test]"
```

## Basic Usage

```python
import numpy as np
from rd2d import rd2d

rng = np.random.default_rng(123)
n = 800
x1 = rng.normal(size=n)
x2 = rng.normal(size=n)
assignment = (x1 >= 0).astype(float)
y = 3 + 2 * x1 + 1.5 * x2 + assignment + rng.normal(size=n)
x = np.column_stack([x1, x2])
b = np.array([[0.0, 0.0], [0.0, 1.0]])

fit = rd2d(y, x, assignment, b, h=0.9, params_cov="main")
fit.main
fit.summary(cbands="main", repp=999).tables["main"]
```

For distance-based designs, pass one signed-distance column per evaluation
point. Nonnegative distances identify observations on the treated side.

```python
from rd2d import rd2d_dist

distance = x1.reshape(-1, 1)
fit_dist = rd2d_dist(y, distance, h=0.5, b=np.array([[0.0, 0.0]]))
fit_dist.main
```

## Development

From this directory:

```sh
python -m pytest
```

## Publishing From GitHub

The workflow `.github/workflows/python-publish.yml` builds, checks, tests, and
publishes the Python package to PyPI with trusted publishing. Configure:

- PyPI project: `rd2d`.
- Trusted publisher owner: `rdpackages`.
- Trusted publisher repository: `rd2d`.
- Trusted publisher workflow: `python-publish.yml`.
- Trusted publisher environment: `pypi`.

In GitHub, create an environment named `pypi`. Add required reviewers there if
you want each PyPI upload to require manual approval. Publish by creating a
GitHub Release or by running the `Publish Python package` workflow manually.
