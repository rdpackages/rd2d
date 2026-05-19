import numpy as np
import pandas as pd

from rd2d import (
    rdbw2d,
    rdbw2d_dist,
    rdbw2d_distance,
    rd2d,
    rd2d_dist,
    rd2d_distance,
    summary,
)


def make_location_data(n=350, seed=101):
    rng = np.random.default_rng(seed)
    x1 = rng.normal(size=n)
    x2 = rng.normal(size=n)
    assignment = (x1 >= 0).astype(float)
    fuzzy = (rng.uniform(size=n) < np.where(assignment == 1, 0.82, 0.18)).astype(float)
    y = 1 + x1 + 0.5 * x2 + assignment + rng.normal(scale=0.5, size=n)
    y_fuzzy = 1 + x1 + 0.5 * x2 + 1.4 * fuzzy + rng.normal(scale=0.5, size=n)
    b = np.array([[0.0, -0.4], [0.0, 0.0], [0.0, 0.4]])
    return y, y_fuzzy, np.column_stack([x1, x2]), assignment, fuzzy, b


def make_distance_data(n=350, seed=202):
    rng = np.random.default_rng(seed)
    x = rng.uniform(-1, 1, size=n)
    assignment = (x >= 0).astype(float)
    fuzzy = (rng.uniform(size=n) < np.where(assignment == 1, 0.8, 0.2)).astype(float)
    y = 1 + 2 * x + 2.5 * assignment + rng.normal(scale=0.3, size=n)
    y_fuzzy = 1 + 2 * x + 1.5 * fuzzy + rng.normal(scale=0.3, size=n)
    return y, y_fuzzy, x.reshape(-1, 1), fuzzy


def test_location_exports_and_tables():
    y, _, x, z, _, b = make_location_data()
    bw = rdbw2d(y, x, z, b, masspoints="off")
    assert list(bw.bws.columns) == ["b1", "b2", "h01", "h02", "h11", "h12", "N.Co", "N.Tr"]

    fit = rd2d(y, x, z, b, h=0.75, masspoints="off", params_cov="main", repp=49)
    assert list(fit.main.columns) == [
        "b1",
        "b2",
        "estimate.p",
        "std.err.p",
        "estimate.q",
        "std.err.q",
        "t.value",
        "p.value",
        "ci.lower",
        "ci.upper",
        "h01",
        "h02",
        "h11",
        "h12",
        "N.Co",
        "N.Tr",
    ]
    assert fit.params_cov["main"].shape == (len(b), len(b))
    summ = summary(fit, cbands="main", WBATE=np.ones(len(b)), LBATE=True)
    assert "main" in summ.tables
    assert list(summ.tables["main"].tail(2).index) == ["WBATE", "LBATE"]


def test_location_fuzzy_tables():
    _, y_fuzzy, x, z, fuzzy, b = make_location_data(seed=303)
    fit = rd2d(
        y_fuzzy,
        x,
        z,
        b,
        h=0.8,
        fuzzy=fuzzy,
        params_other=["itt.0", "fs.1"],
        params_cov=["main", "itt", "fs"],
        masspoints="off",
    )
    assert isinstance(fit.main, pd.DataFrame)
    assert isinstance(fit.itt, pd.DataFrame)
    assert isinstance(fit.fs, pd.DataFrame)
    assert set(fit.params_cov) == {"main", "itt", "fs"}


def test_distance_exports_and_tables():
    y, _, distance, _ = make_distance_data()
    b = np.array([[0.0, 0.0]])
    assert rdbw2d_dist is not None
    assert rdbw2d_distance is not None
    bw = rdbw2d_distance(y, distance, b=b, masspoints="off")
    assert list(bw.bws.columns) == ["b1", "b2", "h0", "h1", "N.Co", "N.Tr"]

    fit = rd2d_distance(y, distance, h=0.45, b=b, p=1, q=1, cbands=False, params_cov="main")
    assert rd2d_dist is not None
    assert list(fit.main.columns) == [
        "b1",
        "b2",
        "estimate.p",
        "std.err.p",
        "estimate.q",
        "std.err.q",
        "t.value",
        "p.value",
        "ci.lower",
        "ci.upper",
        "h0",
        "h1",
        "h0.rbc",
        "h1.rbc",
        "N.Co",
        "N.Tr",
    ]
    assert fit.params_cov["main"].shape == (1, 1)


def test_distance_fuzzy_summary():
    _, y_fuzzy, distance, fuzzy = make_distance_data(seed=404)
    b = np.array([[0.0, 0.0]])
    fit = rd2d_distance(
        y_fuzzy,
        distance,
        h=0.5,
        b=b,
        fuzzy=fuzzy,
        cbands=True,
        params_cov=["itt", "fs"],
        masspoints="off",
        repp=49,
    )
    assert isinstance(fit.itt, pd.DataFrame)
    assert isinstance(fit.fs, pd.DataFrame)
    summ = fit.summary(cbands="main", WBATE=[1.0], LBATE=True)
    assert list(summ.tables["main"].tail(2).index) == ["WBATE", "LBATE"]
