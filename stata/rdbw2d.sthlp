{smcl}
{* *!version 1.0.1  2026-07-08}{...}
{viewerjumpto "Syntax" "rdbw2d##syntax"}{...}
{viewerjumpto "Description" "rdbw2d##description"}{...}
{viewerjumpto "Options" "rdbw2d##options"}{...}
{viewerjumpto "Examples" "rdbw2d##examples"}{...}
{viewerjumpto "Stored results" "rdbw2d##stored_results"}{...}

{title:Title}

{p 4 8}{cmd:rdbw2d} {hline 2} Bandwidth selection for location-based boundary discontinuity designs.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdbw2d} {it:y x1 x2 assignment} {ifin}{cmd:,}
{cmd:b(}{it:# # [# # ...]}{cmd:)}
[{cmd:fuzzy(}{it:varname}{cmd:)}
{cmd:deriv(}{it:# [#]}{cmd:)}
{cmd:tangvec(}{it:# # [# # ...]}{cmd:)}
{cmd:p(}{it:#}{cmd:)}
{cmd:kernel(}{it:kernel}{cmd:)}
{cmd:kerneltype(}{it:type}{cmd:)}
{cmd:bwselect(}{it:selector}{cmd:)}
{cmd:bwparam(}{it:target}{cmd:)}
{cmd:method(}{it:method}{cmd:)}
{cmd:vce(}{it:vcetype}{cmd:)}
{cmd:fitmethod(}{it:method}{cmd:)}
{cmd:covseff(}{it:varlist}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:masspoints(}{it:masspointsoption}{cmd:)}
{cmd:cluster(}{it:clustvar}{cmd:)}
{cmd:scaleregul(}{it:#}{cmd:)}
{cmd:scalebiascrct(}{it:#}{cmd:)}
{cmd:stdvars(}{it:on|off}{cmd:)}]{p_end}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rdbw2d} computes bandwidths for {help rd2d:rd2d}. The returned table has one row per boundary point and columns {cmd:b1}, {cmd:b2}, {cmd:h01}, {cmd:h02}, {cmd:h11}, {cmd:h12}, {cmd:N_Co}, and {cmd:N_Tr}.{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{cmd:b()} specifies boundary evaluation points as pairs and is required.{p_end}

{p 4 8}{cmd:fuzzy()} specifies treatment receipt/status for fuzzy bandwidth selection. If supplied, {cmd:bwparam()} controls whether the selector targets the fuzzy Wald ratio or the reduced-form outcome discontinuity.{p_end}

{p 4 8}{cmd:deriv()} specifies derivative orders in the first and second coordinates. Default is {cmd:deriv(0 0)}.{p_end}

{p 4 8}{cmd:tangvec()} supplies directional derivative vectors, one pair per boundary point. If supplied, {cmd:tangvec()} overrides {cmd:deriv()} and requires {cmd:p(1)} or higher.{p_end}

{p 4 8}{cmd:p()} sets the local polynomial order. Default is {cmd:p(1)}.{p_end}

{p 4 8}{cmd:kernel()} specifies the kernel function. Available choices are:{p_end}
{p 8 12}{cmd:triangular} or {cmd:tri}: triangular kernel, the default.{p_end}
{p 8 12}{cmd:epanechnikov} or {cmd:epa}: Epanechnikov kernel.{p_end}
{p 8 12}{cmd:uniform} or {cmd:uni}: uniform kernel.{p_end}
{p 8 12}{cmd:gaussian} or {cmd:gau}: Gaussian kernel.{p_end}

{p 4 8}{cmd:kerneltype()} specifies the kernel structure. Available choices are:{p_end}
{p 8 12}{cmd:prod}: product kernel, the default.{p_end}
{p 8 12}{cmd:rad}: radial kernel.{p_end}

{p 4 8}{cmd:bwselect()} specifies the bandwidth selector. Available choices are:{p_end}
{p 8 12}{cmd:mserd}: one common MSE-optimal bandwidth selector for the boundary RD treatment-effect estimator at each evaluation point, the default.{p_end}
{p 8 12}{cmd:cerrd}: CER-optimal counterpart of {cmd:mserd}.{p_end}
{p 8 12}{cmd:imserd}: IMSE-optimal common bandwidth selector for the boundary RD treatment-effect estimator based on all evaluation points.{p_end}
{p 8 12}{cmd:icerrd}: CER-optimal counterpart of {cmd:imserd}.{p_end}
{p 8 12}{cmd:msetwo}: two different MSE-optimal bandwidth selectors, one for the control group and one for the treated group, at each evaluation point.{p_end}
{p 8 12}{cmd:certwo}: CER-optimal counterpart of {cmd:msetwo}.{p_end}
{p 8 12}{cmd:imsetwo}: two IMSE-optimal bandwidth selectors, one for each group, based on all evaluation points.{p_end}
{p 8 12}{cmd:icertwo}: CER-optimal counterpart of {cmd:imsetwo}.{p_end}

{p 4 8}{cmd:bwparam()} specifies the target parameter for fuzzy bandwidth selection. Available choices are:{p_end}
{p 8 12}{cmd:main}: selects bandwidths for the linearized fuzzy Wald ratio, the default.{p_end}
{p 8 12}{cmd:itt}: selects bandwidths for the reduced-form outcome discontinuity. This option is ignored in sharp designs.{p_end}

{p 4 8}{cmd:method()} specifies the bandwidth calculation method. Available choices are:{p_end}
{p 8 12}{cmd:dpi}: data-driven plug-in MSE-optimal selector, the default.{p_end}
{p 8 12}{cmd:rot}: rule-of-thumb selector.{p_end}

{p 4 8}{cmd:vce()} specifies the variance-covariance estimator used in the bandwidth calculations. Available choices are {cmd:hc0}, {cmd:hc1}, {cmd:hc2}, and {cmd:hc3}. Default is {cmd:hc1}.{p_end}

{p 4 8}{cmd:fitmethod()} specifies the fitting convention used for variance constants in bandwidth selection. The default, {cmd:fitmethod(joint)}, uses the joint treatment-interacted regression convention. {cmd:fitmethod(separate)} uses the legacy two-sample side-specific convention.{p_end}

{p 4 8}{cmd:covseff(}{it:varlist}{cmd:)} specifies pre-intervention covariates used for efficiency adjustment in bandwidth selection. Covariates enter with common coefficients across treatment sides after residualizing on the side-specific local polynomial bases.{p_end}

{p 4 8}{cmd:bwcheck()} enlarges preliminary bandwidths, when needed, so that at least the requested number of unique observations are included in each side-specific window. Default is {cmd:bwcheck(20)}.{p_end}

{p 4 8}{cmd:masspoints()} specifies how repeated running-variable values are handled. Available choices are:{p_end}
{p 8 12}{cmd:check}: detects mass points and reports unique-observation counts, the default.{p_end}
{p 8 12}{cmd:adjust}: adjusts preliminary bandwidths to ensure a minimum number of unique observations within each side of the boundary.{p_end}
{p 8 12}{cmd:off}: ignores mass points.{p_end}

{p 4 8}{cmd:cluster()} specifies a cluster ID variable for cluster-robust bandwidth calculations.{p_end}

{p 4 8}{cmd:scaleregul()} specifies the scaling factor for the regularization term in bandwidth selection. Default is {cmd:scaleregul(1)}.{p_end}

{p 4 8}{cmd:scalebiascrct()} specifies the scaling factor used for bias correction based on higher-order expansions. Default is {cmd:scalebiascrct(1)}.{p_end}

{p 4 8}{cmd:stdvars()} controls whether the running variables are standardized before computing automatic bandwidths. Use {cmd:stdvars(on)} or {cmd:stdvars(off)}. Default is {cmd:on}.{p_end}

{marker examples}{...}
{title:Examples}

{p 8 8}{cmd:. clear}{p_end}
{p 8 8}{cmd:. set type double}{p_end}
{p 8 8}{cmd:. set obs 800}{p_end}
{p 8 8}{cmd:. set seed 123}{p_end}
{p 8 8}{cmd:. generate double x1 = rnormal()}{p_end}
{p 8 8}{cmd:. generate double x2 = rnormal()}{p_end}
{p 8 8}{cmd:. generate double d = x1 >= 0}{p_end}
{p 8 8}{cmd:. generate double y = 3 + 2*x1 + 1.5*x2 + d + rnormal()}{p_end}
{p 8 8}{cmd:. rdbw2d y x1 x2 d, b(0 0 0 1)}{p_end}
{p 8 8}{cmd:. generate double z = x1 + x2 + rnormal()}{p_end}
{p 8 8}{cmd:. rdbw2d y x1 x2 d, b(0 0 0 1) covseff(z) fitmethod(joint)}{p_end}

{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:rdbw2d} stores {cmd:e(bws)}, {cmd:e(N)}, {cmd:e(p)}, and option labels in {cmd:e()}, including {cmd:e(fitmethod)} and {cmd:e(covseff)} when supplied.{p_end}
