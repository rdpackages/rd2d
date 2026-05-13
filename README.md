# Boundary Discontinuity Designs

[![R-CMD-check](https://github.com/rdpackages/rd2d/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rdpackages/rd2d/actions/workflows/R-CMD-check.yaml)

The `rd2d` package provides R implementations of pointwise and uniform estimation and inference for Boundary Discontinuity Designs employing local polynomial methods.
The main functions are `rd2d()` and `rdbw2d()` for location-based methods, and `rd2d.distance()` and `rdbw2d.distance()` for distance-based methods. Both approaches support sharp and fuzzy designs. The distance-based methods target the level of the boundary average treatment effect curve using bivariate scores reduced to signed distance-based running variables.
Summary methods can also report weighted and largest boundary average treatment effects (WBATE and LBATE) when the fitted object stores the required covariance matrix, including intention-to-treat, first-stage, and fuzzy outputs.

This work was supported in part by the National Science Foundation through grants [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432), [DMS-2210561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2210561), [SES-2241575](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2241575), and [SES-2342226](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2342226), and by the National Institute for Food and Agriculture through grant [2024-67023-42704](https://www.nifa.usda.gov/data).


## R Implementation

To install/update in R type:
```
install.packages('rd2d')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rd2d/rd2d.pdf), [CRAN repository](https://cran.r-project.org/package=rd2d).

- Illustration: [Simulation and estimation](R/rd2d_illustration.R), [Plots](R/rd2d_plot.R). Generated illustration output is written to `R/output/` and is not tracked by git.

## Development

The R package source lives in [`R/rd2d`](R/rd2d/). The remaining files in [`R`](R/) are illustration scripts and generated output.

From the repository root, check the package with:

```sh
R CMD check --no-manual R/rd2d
```

GitHub Actions runs this check on Linux, macOS, and Windows for pull requests and pushes that touch the R package source.

Research papers and replication files are maintained separately from this package repository.


## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Cattaneo, Titiunik and Yu (2025): [rd2d: Causal Inference in Boundary Discontinuity Designs](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2025_rd2d.pdf).<br>
Working paper.


### Technical and Methodological

- Cattaneo, Titiunik and Yu (2026): [Boundary Discontinuity Designs: Theory and Practice](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_ESWC.pdf).<br>
Advances in Economics and Econometrics: Thirteenth World Congress, eds. R. Griffith, Y. Gorodnichenko, M. Kandori and F. Molinari, Cambridge University Press, Vol. 1, Ch. 2, to appear, 2026.

- Cattaneo, Titiunik and Yu (2026): [Estimation and Inference in Boundary Discontinuity Designs: Location-Based Methods](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_JASA.pdf).<br>
_Journal of the American Statistical Association_, revise and resubmit.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_JASA--Supplement.pdf).

- Cattaneo, Titiunik and Yu (2026): [Estimation and Inference in Boundary Discontinuity Designs: Distance-Based Methods](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_JOE.pdf).<br>
_Journal of Econometrics_, revise and resubmit.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_JOE--Supplement.pdf).

- Cattaneo, Titiunik and Yu (2026): [Estimation and Inference in Boundary Discontinuity Designs: Pooling-Based Methods](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_BDD-Pooling.pdf).<br>
Working paper.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Cattaneo-Titiunik-Yu_2026_BDD-Pooling--Supplement.pdf).


<br><br>
