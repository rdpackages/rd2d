import numpy as np
import pandas as pd
import pytest

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


def make_covariates(x, seed=505):
    rng = np.random.default_rng(seed)
    z1 = 0.4 * x[:, 0] - 0.2 * x[:, 1] + rng.normal(scale=0.4, size=x.shape[0])
    z2 = rng.normal(size=x.shape[0])
    return np.column_stack([z1, z2])


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

    fit = rd2d(y, x, z, b, h=0.75, masspoints="off", params_cov="main")
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
    assert "ci" in fit
    assert "cb" not in fit
    summ = summary(fit, cbands="main", repp=49, WBATE=np.ones(len(b)), LBATE=True)
    assert "main" in summ.tables
    assert {"ci.lower", "ci.upper", "cb.lower", "cb.upper"}.issubset(summ.tables["main"].columns)
    assert list(summ.tables["main"].tail(2).index) == ["WBATE", "LBATE"]
    with np.testing.assert_raises_regex(ValueError, "repp must be a positive integer"):
        summary(fit, cbands="main", repp=0)


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
    )
    assert isinstance(fit.itt, pd.DataFrame)
    assert isinstance(fit.fs, pd.DataFrame)
    summ = fit.summary(cbands="main", repp=49, WBATE=[1.0], LBATE=True)
    assert "ci" in fit
    assert "cb" not in fit
    assert {"ci.lower", "ci.upper", "cb.lower", "cb.upper"}.issubset(summ.tables["main"].columns)
    assert list(summ.tables["main"].tail(2).index) == ["WBATE", "LBATE"]


def test_location_fitmethod_and_covariates():
    y, _, x, z, _, b = make_location_data(seed=515)
    covs = make_covariates(x)
    fit_joint = rd2d(y, x, z, b, h=0.75, vce="hc0", masspoints="off", fitmethod="joint")
    fit_sep = rd2d(y, x, z, b, h=0.75, vce="hc0", masspoints="off", fitmethod="separate")
    np.testing.assert_allclose(fit_joint.main["estimate.p"], fit_sep.main["estimate.p"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(fit_joint.main["std.err.q"], fit_sep.main["std.err.q"], rtol=1e-12, atol=1e-12)

    fit_cov = rd2d(y, x, z, b, h=0.75, covs_eff=covs, masspoints="off")
    assert fit_cov.opt["covs.eff"] is True
    assert fit_cov.opt["N.covs.eff"] == 2
    assert fit_cov.opt["fitmethod"] == "joint"
    assert np.all(np.isfinite(fit_cov.main["estimate.p"]))

    with pytest.warns(RuntimeWarning, match="rank deficient"):
        fit_rank = rd2d(y, x, z, b, h=0.75, covs_eff=np.column_stack([covs[:, 0], covs[:, 0]]), masspoints="off")
    assert fit_rank.opt["covs.rank.deficient"] is True
    assert fit_rank.opt["N.covs.used"] == 1


def test_location_cluster_separate_covariance_is_side_additive():
    y, _, x, z, _, b = make_location_data(n=420, seed=818)
    cluster = np.arange(len(y)) % 35
    fit = rd2d(
        y,
        x,
        z,
        b,
        h=0.8,
        vce="hc1",
        cluster=cluster,
        fitmethod="separate",
        masspoints="off",
        params_cov=["main", "main.0", "main.1"],
    )
    np.testing.assert_allclose(
        fit.params_cov["main"],
        fit.params_cov["main.0"] + fit.params_cov["main.1"],
        rtol=1e-12,
        atol=1e-12,
    )


def test_location_cluster_labels_are_partition_invariant():
    y, _, x, z, _, b = make_location_data(n=420, seed=828)
    codes = (np.arange(len(y)) * 7 + 11) % 43
    numeric_cluster = codes * 10 + 3
    string_cluster = np.asarray([f"group-{1000 - code}" for code in codes], dtype=object)

    for fitmethod in ["separate", "joint"]:
        fit_numeric = rd2d(
            y,
            x,
            z,
            b,
            h=0.8,
            vce="hc1",
            cluster=numeric_cluster,
            fitmethod=fitmethod,
            masspoints="off",
            params_cov="main",
        )
        fit_string = rd2d(
            y,
            x,
            z,
            b,
            h=0.8,
            vce="hc1",
            cluster=string_cluster,
            fitmethod=fitmethod,
            masspoints="off",
            params_cov="main",
        )
        np.testing.assert_allclose(fit_numeric.main.to_numpy(float), fit_string.main.to_numpy(float), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fit_numeric.params_cov["main"], fit_string.params_cov["main"], rtol=1e-12, atol=1e-12)


def test_distance_fitmethod_and_covariates():
    y, _, distance, _ = make_distance_data(seed=616)
    b = np.array([[0.0, 0.0]])
    rng = np.random.default_rng(717)
    covs = np.column_stack([rng.normal(size=len(y)), rng.normal(size=len(y))])
    fit_joint = rd2d_distance(y, distance, h=0.5, b=b, vce="hc0", fitmethod="joint", cbands=False)
    fit_sep = rd2d_distance(y, distance, h=0.5, b=b, vce="hc0", fitmethod="separate", cbands=False)
    np.testing.assert_allclose(fit_joint.main["estimate.q"], fit_sep.main["estimate.q"], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(fit_joint.main["std.err.q"], fit_sep.main["std.err.q"], rtol=1e-12, atol=1e-12)

    fit_cov = rd2d_distance(y, distance, h=0.5, b=b, covs_eff=covs, cbands=False)
    assert fit_cov.opt["covs.eff"] is True
    assert fit_cov.opt["N.covs.eff"] == 2
    assert np.all(np.isfinite(fit_cov.main["estimate.p"]))

    bw_cov = rdbw2d_distance(y, distance, b=b, covs_eff=covs, masspoints="off")
    assert bw_cov.opt["covs.eff"] is True
    assert list(bw_cov.bws.columns) == ["b1", "b2", "h0", "h1", "N.Co", "N.Tr"]


def test_distance_cluster_separate_covariance_is_side_additive():
    y, _, distance, _ = make_distance_data(n=420, seed=919)
    b = np.array([[0.0, 0.0]])
    cluster = np.arange(len(y)) % 35
    fit = rd2d_distance(
        y,
        distance,
        h=0.5,
        b=b,
        vce="hc1",
        cluster=cluster,
        fitmethod="separate",
        cbands=False,
        params_cov=["main", "main.0", "main.1"],
    )
    np.testing.assert_allclose(
        fit.params_cov["main"],
        fit.params_cov["main.0"] + fit.params_cov["main.1"],
        rtol=1e-12,
        atol=1e-12,
    )


def test_distance_cluster_labels_are_partition_invariant():
    y, _, distance, _ = make_distance_data(n=420, seed=929)
    b = np.array([[0.0, 0.0]])
    codes = (np.arange(len(y)) * 11 + 5) % 37
    numeric_cluster = codes * 100 + 17
    string_cluster = np.asarray([f"site-{500 - code}" for code in codes], dtype=object)

    for fitmethod in ["separate", "joint"]:
        fit_numeric = rd2d_distance(
            y,
            distance,
            h=0.5,
            b=b,
            vce="hc1",
            cluster=numeric_cluster,
            fitmethod=fitmethod,
            cbands=False,
            params_cov="main",
        )
        fit_string = rd2d_distance(
            y,
            distance,
            h=0.5,
            b=b,
            vce="hc1",
            cluster=string_cluster,
            fitmethod=fitmethod,
            cbands=False,
            params_cov="main",
        )
        np.testing.assert_allclose(fit_numeric.main.to_numpy(float), fit_string.main.to_numpy(float), rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(fit_numeric.params_cov["main"], fit_string.params_cov["main"], rtol=1e-12, atol=1e-12)
