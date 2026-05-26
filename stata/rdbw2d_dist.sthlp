{smcl}
{* *!version 1.0.0  2026-05-26}{...}
{viewerjumpto "Syntax" "rdbw2d_dist##syntax"}{...}
{viewerjumpto "Description" "rdbw2d_dist##description"}{...}
{viewerjumpto "Options" "rdbw2d_dist##options"}{...}
{viewerjumpto "Examples" "rdbw2d_dist##examples"}{...}
{viewerjumpto "Stored results" "rdbw2d_dist##stored_results"}{...}

{title:Title}

{p 4 8}{cmd:rdbw2d_dist} {hline 2} Bandwidth selection for distance-based boundary discontinuity designs.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdbw2d_dist} {it:y distancevar [distancevar ...]} {ifin}
[{cmd:,} {cmd:b(}{it:# # [# # ...]}{cmd:)}
{cmd:fuzzy(}{it:varname}{cmd:)}
{cmd:p(}{it:#}{cmd:)}
{cmd:kinkunknown(}{it:on|off|#}{cmd:)}
{cmd:kinkposition(}{it:# [# ...]}{cmd:)}
{cmd:kernel(}{it:kernel}{cmd:)}
{cmd:bwselect(}{it:selector}{cmd:)}
{cmd:bwparam(}{it:target}{cmd:)}
{cmd:vce(}{it:vcetype}{cmd:)}
{cmd:fitmethod(}{it:method}{cmd:)}
{cmd:covseff(}{it:varlist}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:masspoints(}{it:masspointsoption}{cmd:)}
{cmd:cluster(}{it:clustvar}{cmd:)}
{cmd:scaleregul(}{it:#}{cmd:)}
{cmd:cqt(}{it:#}{cmd:)}]{p_end}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rdbw2d_dist} computes distance-based bandwidths for {help rd2d_dist:rd2d_dist}. The alias {cmd:rdbw2d_distance} is also provided.{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{cmd:b()} supplies optional boundary coordinates, one pair per distance variable.{p_end}

{p 4 8}{cmd:fuzzy()} specifies treatment receipt/status for fuzzy bandwidth selection. If supplied, {cmd:bwparam()} controls whether the selector targets the fuzzy Wald ratio or the reduced-form outcome discontinuity.{p_end}

{p 4 8}{cmd:p()} sets the local polynomial order. Default is {cmd:p(1)}.{p_end}

{p 4 8}{cmd:kinkunknown()} requests the unknown-kink bandwidth-rate adjustment. Use {cmd:on}, {cmd:off}, or a 0/1 indicator. For bandwidth selection, the first unknown-kink indicator is used.{p_end}

{p 4 8}{cmd:kinkposition()} supplies known kink locations for adaptive automatic bandwidth selection. Provide 1-based boundary-point indices, or a 0/1 indicator with one value per distance variable. This option requires {cmd:b()} and cannot be combined with {cmd:kinkunknown()}.{p_end}

{p 4 8}{cmd:kernel()} specifies the kernel function. Available choices are:{p_end}
{p 8 12}{cmd:triangular} or {cmd:tri}: triangular kernel, the default.{p_end}
{p 8 12}{cmd:epanechnikov} or {cmd:epa}: Epanechnikov kernel.{p_end}
{p 8 12}{cmd:uniform} or {cmd:uni}: uniform kernel.{p_end}
{p 8 12}{cmd:gaussian} or {cmd:gau}: Gaussian kernel.{p_end}

{p 4 8}{cmd:bwselect()} specifies the bandwidth selector. Available choices are:{p_end}
{p 8 12}{cmd:mserd}: one common MSE-optimal bandwidth selector for the boundary RD treatment-effect estimator at each evaluation point, the default.{p_end}
{p 8 12}{cmd:cerrd}: CER-optimal counterpart of {cmd:mserd}.{p_end}
{p 8 12}{cmd:imserd}: IMSE-optimal common bandwidth selector for the boundary RD treatment-effect estimator based on all evaluation points.{p_end}
{p 8 12}{cmd:icerrd}: CER-optimal counterpart of {cmd:imserd}.{p_end}
{p 8 12}{cmd:msetwo}: two different MSE-optimal bandwidth selectors, one for negative signed distances and one for nonnegative signed distances, at each evaluation point.{p_end}
{p 8 12}{cmd:certwo}: CER-optimal counterpart of {cmd:msetwo}.{p_end}
{p 8 12}{cmd:imsetwo}: two IMSE-optimal bandwidth selectors, one for each side, based on all evaluation points.{p_end}
{p 8 12}{cmd:icertwo}: CER-optimal counterpart of {cmd:imsetwo}.{p_end}

{p 4 8}{cmd:bwparam()} specifies the target parameter for fuzzy bandwidth selection. Available choices are:{p_end}
{p 8 12}{cmd:main}: selects bandwidths for the linearized fuzzy Wald ratio, the default.{p_end}
{p 8 12}{cmd:itt}: selects bandwidths for the reduced-form outcome discontinuity. This option is ignored in sharp designs.{p_end}

{p 4 8}{cmd:vce()} specifies the variance-covariance estimator used in bandwidth calculations. Available choices are {cmd:hc0}, {cmd:hc1}, {cmd:hc2}, and {cmd:hc3}. Default is {cmd:hc1}.{p_end}

{p 4 8}{cmd:fitmethod()} specifies the fitting convention used for variance constants in bandwidth selection. The default, {cmd:fitmethod(joint)}, uses the joint treatment-interacted regression convention. {cmd:fitmethod(separate)} uses the legacy two-sample side-specific convention.{p_end}

{p 4 8}{cmd:covseff(}{it:varlist}{cmd:)} specifies pre-intervention covariates used for efficiency adjustment in bandwidth selection. Covariates enter with common coefficients across treatment sides after residualizing on the side-specific local polynomial bases.{p_end}

{p 4 8}{cmd:bwcheck()} enlarges preliminary bandwidths, when needed, so that at least the requested number of unique observations are included in each side-specific window. Default is {cmd:20 + p + 1}.{p_end}

{p 4 8}{cmd:masspoints()} specifies how repeated signed-distance values are handled. Available choices are:{p_end}
{p 8 12}{cmd:check}: detects mass points and reports unique-observation counts, the default.{p_end}
{p 8 12}{cmd:adjust}: adjusts preliminary bandwidths to ensure a minimum number of unique observations within each side of the cutoff.{p_end}
{p 8 12}{cmd:off}: ignores mass points.{p_end}

{p 4 8}{cmd:cluster()} specifies a cluster ID variable for cluster-robust bandwidth calculations.{p_end}

{p 4 8}{cmd:scaleregul()} specifies the scaling factor for the regularization term in bandwidth selection. Default is {cmd:scaleregul(1)}.{p_end}

{p 4 8}{cmd:cqt()} specifies the subsample quantile used for initial bias estimation in automatic bandwidth selection. Default is {cmd:cqt(.5)}.{p_end}

{marker examples}{...}
{title:Examples}

{p 8 8}{cmd:. clear}{p_end}
{p 8 8}{cmd:. set type double}{p_end}
{p 8 8}{cmd:. set obs 800}{p_end}
{p 8 8}{cmd:. set seed 123}{p_end}
{p 8 8}{cmd:. generate double x = runiform()*2 - 1}{p_end}
{p 8 8}{cmd:. generate double d = x >= 0}{p_end}
{p 8 8}{cmd:. generate double y = 1 + 2*x + 1.5*d + rnormal()*.5}{p_end}
{p 8 8}{cmd:. rdbw2d_dist y x, b(0 0)}{p_end}
{p 8 8}{cmd:. generate double z = x + rnormal()}{p_end}
{p 8 8}{cmd:. rdbw2d_dist y x, b(0 0) covseff(z) fitmethod(joint)}{p_end}

{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:rdbw2d_dist} stores {cmd:e(bws)}, {cmd:e(N)}, {cmd:e(p)}, and option labels in {cmd:e()}, including {cmd:e(fitmethod)} and {cmd:e(covseff)} when supplied.{p_end}
