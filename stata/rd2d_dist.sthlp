{smcl}
{* *!version 1.0.0  2026-05-26}{...}
{viewerjumpto "Syntax" "rd2d_dist##syntax"}{...}
{viewerjumpto "Description" "rd2d_dist##description"}{...}
{viewerjumpto "Options" "rd2d_dist##options"}{...}
{viewerjumpto "Examples" "rd2d_dist##examples"}{...}
{viewerjumpto "Stored results" "rd2d_dist##stored_results"}{...}

{title:Title}

{p 4 8}{cmd:rd2d_dist} {hline 2} Distance-based methods for boundary discontinuity designs.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rd2d_dist} {it:y distancevar [distancevar ...]} {ifin}
[{cmd:,} {cmd:b(}{it:# # [# # ...]}{cmd:)}
{cmd:h(}{it:# | # # [# # ...]}{cmd:)}
{cmd:fuzzy(}{it:varname}{cmd:)}
{cmd:p(}{it:#}{cmd:)}
{cmd:q(}{it:#}{cmd:)}
{cmd:kinkunknown(}{it:on|off|# #}{cmd:)}
{cmd:kinkposition(}{it:# [# ...]}{cmd:)}
{cmd:kernel(}{it:kernel}{cmd:)}
{cmd:level(}{it:#}{cmd:)}
{cmd:side(}{it:side}{cmd:)}
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

{p 4 8}{cmd:rd2d_dist} estimates distance-based local polynomial boundary discontinuity effects using signed distance scores. Negative distances identify the control side and nonnegative distances identify the treated side. Each distance variable corresponds to one boundary evaluation point.{p_end}

{p 4 8}The alias {cmd:rd2d_distance} is also provided.{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Boundary Points and Design}

{p 4 8}{cmd:b()} optionally supplies two-dimensional boundary point coordinates, one pair per distance variable. If omitted, coordinates are stored as missing in result tables.{p_end}

{p 4 8}{cmd:h()} supplies user-provided bandwidths. A scalar uses the same bandwidth for both groups and all evaluation points. Otherwise, supply two values for each evaluation point: control bandwidth {it:h0} and treated bandwidth {it:h1}. If {cmd:h()} is omitted, {cmd:rd2d_dist} computes automatic bandwidths using {help rdbw2d_dist:rdbw2d_dist}.{p_end}

{p 4 8}{cmd:fuzzy()} specifies treatment receipt/status for fuzzy designs. The main table reports the fuzzy Wald ratio. The reduced-form outcome discontinuity is stored in {cmd:e(itt)}, and the first-stage treatment receipt/status discontinuity is stored in {cmd:e(fs)}. Side-specific fuzzy tables are stored in {cmd:e(itt_0)}, {cmd:e(itt_1)}, {cmd:e(fs_0)}, and {cmd:e(fs_1)}.{p_end}

{dlgtab:Local Polynomial Regression}

{p 4 8}{cmd:p()} and {cmd:q()} set the local polynomial orders for point estimation and bias-corrected inference. Defaults are {cmd:p(1)} and {cmd:q(p+1)} unless unknown-kink inference is requested, in which case the default inference order follows the unknown-kink setting.{p_end}

{p 4 8}{cmd:kinkunknown()} controls unknown-kink bandwidth-rate adjustments. Use {cmd:on} to apply both point-estimation and inference adjustments, {cmd:off} to apply neither adjustment, or two 0/1 indicators where the first controls point-estimation bandwidths and the second controls inference bandwidths. The second indicator can be 1 only when the first is 1.{p_end}

{p 4 8}{cmd:kinkposition()} supplies known kink locations for adaptive automatic bandwidth selection. Provide 1-based boundary-point indices, or a 0/1 indicator with one value per distance variable. This option requires {cmd:b()}, applies only when {cmd:h()} is omitted, and cannot be combined with {cmd:kinkunknown()}.{p_end}

{p 4 8}{cmd:kernel()} specifies the kernel function. Available choices are:{p_end}
{p 8 12}{cmd:triangular} or {cmd:tri}: triangular kernel, the default.{p_end}
{p 8 12}{cmd:epanechnikov} or {cmd:epa}: Epanechnikov kernel.{p_end}
{p 8 12}{cmd:uniform} or {cmd:uni}: uniform kernel.{p_end}
{p 8 12}{cmd:gaussian} or {cmd:gau}: Gaussian kernel.{p_end}

{dlgtab:Inference}

{p 4 8}{cmd:vce()} specifies the variance-covariance estimator. Available choices are:{p_end}
{p 8 12}{cmd:hc0}: heteroskedasticity-robust plug-in residual variance estimator without small-sample adjustment.{p_end}
{p 8 12}{cmd:hc1}: heteroskedasticity-robust plug-in residual variance estimator with HC1 degrees-of-freedom adjustment, the default.{p_end}
{p 8 12}{cmd:hc2}: heteroskedasticity-robust plug-in residual variance estimator with HC2 leverage adjustment.{p_end}
{p 8 12}{cmd:hc3}: heteroskedasticity-robust plug-in residual variance estimator with HC3 leverage adjustment.{p_end}

{p 4 8}{cmd:cluster()} specifies a cluster ID variable. When supplied, cluster-level influence sums with degrees-of-freedom weights are used.{p_end}

{p 4 8}{cmd:fitmethod()} specifies the fitting convention used for estimation, inference, and automatic bandwidth selection. The default, {cmd:fitmethod(joint)}, uses the joint treatment-interacted regression convention. {cmd:fitmethod(separate)} uses the legacy two-sample side-specific fitting convention. Point estimates are unchanged without covariates, but HC1 and clustered-robust standard errors may differ because the joint method counts degrees of freedom and clusters jointly across sides.{p_end}

{p 4 8}{cmd:covseff(}{it:varlist}{cmd:)} specifies pre-intervention covariates used for efficiency adjustment. Covariates enter with common coefficients across treatment sides after residualizing on the side-specific local polynomial bases. When {cmd:h()} is omitted, the same covariate adjustment is propagated to automatic bandwidth selection.{p_end}

{p 4 8}{cmd:level()} sets the confidence level. Default is {cmd:level(95)}.{p_end}

{p 4 8}{cmd:side()} specifies the type of confidence interval. Available choices are:{p_end}
{p 8 12}{cmd:two}: two-sided intervals, the default.{p_end}
{p 8 12}{cmd:left}: left-tail one-sided intervals.{p_end}
{p 8 12}{cmd:right}: right-tail one-sided intervals.{p_end}

{dlgtab:Bandwidth Selection}

{p 4 8}{cmd:bwselect()} specifies the automatic bandwidth selector. Available choices are:{p_end}
{p 8 12}{cmd:mserd}: one common MSE-optimal bandwidth selector for the boundary RD treatment-effect estimator at each evaluation point, the default.{p_end}
{p 8 12}{cmd:cerrd}: CER-optimal counterpart of {cmd:mserd}.{p_end}
{p 8 12}{cmd:imserd}: IMSE-optimal common bandwidth selector for the boundary RD treatment-effect estimator based on all evaluation points.{p_end}
{p 8 12}{cmd:icerrd}: CER-optimal counterpart of {cmd:imserd}.{p_end}
{p 8 12}{cmd:msetwo}: two different MSE-optimal bandwidth selectors, one for negative signed distances and one for nonnegative signed distances, at each evaluation point.{p_end}
{p 8 12}{cmd:certwo}: CER-optimal counterpart of {cmd:msetwo}.{p_end}
{p 8 12}{cmd:imsetwo}: two IMSE-optimal bandwidth selectors, one for each side, based on all evaluation points.{p_end}
{p 8 12}{cmd:icertwo}: CER-optimal counterpart of {cmd:imsetwo}.{p_end}
{p 8 12}{cmd:user provided}: used when {cmd:h()} is supplied.{p_end}

{p 4 8}{cmd:bwparam()} specifies the target parameter for fuzzy automatic bandwidth selection. Available choices are:{p_end}
{p 8 12}{cmd:main}: selects bandwidths for the linearized fuzzy Wald ratio, the default.{p_end}
{p 8 12}{cmd:itt}: selects bandwidths for the reduced-form outcome discontinuity. This option is ignored in sharp designs and when {cmd:h()} is supplied.{p_end}

{p 4 8}{cmd:bwcheck()} enlarges preliminary bandwidths, when needed, so that at least the requested number of unique observations are included in each side-specific window. Default is {cmd:50 + p + 1}.{p_end}

{p 4 8}{cmd:masspoints()} specifies how repeated signed-distance values are handled. Available choices are:{p_end}
{p 8 12}{cmd:check}: detects mass points and reports unique-observation counts, the default.{p_end}
{p 8 12}{cmd:adjust}: adjusts preliminary bandwidths to ensure a minimum number of unique observations within each side of the cutoff.{p_end}
{p 8 12}{cmd:off}: ignores mass points.{p_end}

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
{p 8 8}{cmd:. rd2d_dist y x, b(0 0) h(.5)}{p_end}
{p 8 8}{cmd:. generate double z = x + rnormal()}{p_end}
{p 8 8}{cmd:. rd2d_dist y x, b(0 0) h(.5) covseff(z) fitmethod(joint)}{p_end}

{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:rd2d_dist} stores {cmd:e(main)}, {cmd:e(bw)}, {cmd:e(V_main)}, and, when applicable, {cmd:e(itt)}, {cmd:e(fs)}, {cmd:e(itt_0)}, {cmd:e(itt_1)}, {cmd:e(fs_0)}, {cmd:e(fs_1)}, plus covariance matrices {cmd:e(V_itt)}, {cmd:e(V_fs)}, {cmd:e(V_itt_0)}, {cmd:e(V_itt_1)}, {cmd:e(V_fs_0)}, and {cmd:e(V_fs_1)}. Option labels include {cmd:e(fitmethod)} and {cmd:e(covseff)} when supplied.{p_end}
