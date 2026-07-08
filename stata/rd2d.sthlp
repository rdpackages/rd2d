{smcl}
{* *!version 1.0.1  2026-07-08}{...}
{viewerjumpto "Syntax" "rd2d##syntax"}{...}
{viewerjumpto "Description" "rd2d##description"}{...}
{viewerjumpto "Options" "rd2d##options"}{...}
{viewerjumpto "Examples" "rd2d##examples"}{...}
{viewerjumpto "Stored results" "rd2d##stored_results"}{...}
{viewerjumpto "References" "rd2d##references"}{...}
{viewerjumpto "Authors" "rd2d##authors"}{...}

{title:Title}

{p 4 8}{cmd:rd2d} {hline 2} Location-based methods for boundary discontinuity designs.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rd2d} {it:y x1 x2 assignment} {ifin}{cmd:,}
{cmd:b(}{it:# # [# # ...]}{cmd:)}
[{cmd:h(}{it:# | # # # # [# # # # ...]}{cmd:)}
{cmd:fuzzy(}{it:varname}{cmd:)}
{cmd:deriv(}{it:# [#]}{cmd:)}
{cmd:tangvec(}{it:# # [# # ...]}{cmd:)}
{cmd:p(}{it:#}{cmd:)}
{cmd:q(}{it:#}{cmd:)}
{cmd:kernel(}{it:kernel}{cmd:)}
{cmd:kerneltype(}{it:type}{cmd:)}
{cmd:vce(}{it:vcetype}{cmd:)}
{cmd:fitmethod(}{it:method}{cmd:)}
{cmd:covseff(}{it:varlist}{cmd:)}
{cmd:cluster(}{it:clustvar}{cmd:)}
{cmd:level(}{it:#}{cmd:)}
{cmd:side(}{it:side}{cmd:)}
{cmd:bwselect(}{it:selector}{cmd:)}
{cmd:bwparam(}{it:target}{cmd:)}
{cmd:method(}{it:method}{cmd:)}
{cmd:bwcheck(}{it:#}{cmd:)}
{cmd:masspoints(}{it:masspointsoption}{cmd:)}
{cmd:scaleregul(}{it:#}{cmd:)}
{cmd:scalebiascrct(}{it:#}{cmd:)}
{cmd:stdvars(}{it:on|off}{cmd:)}]{p_end}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rd2d} estimates location-based local polynomial boundary discontinuity effects at one or more two-dimensional boundary evaluation points. The first three variables are the outcome and bivariate running variable; {it:assignment} is a binary treatment-side indicator.{p_end}

{p 4 8}The option {cmd:b()} supplies boundary points as pairs. For example, {cmd:b(0 0 0 1)} requests estimates at (0,0) and (0,1). The option {cmd:h()} accepts either a scalar bandwidth or four bandwidths per boundary point: control x1, control x2, treated x1, and treated x2.{p_end}

{p 4 8}Companion commands are {help rdbw2d:rdbw2d} for location-based bandwidth selection and {help rd2d_dist:rd2d_dist} for distance-based designs. Requires Stata 16 or later.{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Boundary Points and Design}

{p 4 8}{cmd:b(}{it:# # [# # ...]}{cmd:)} specifies the boundary evaluation points and is required.{p_end}

{p 4 8}{cmd:h()} specifies user-provided bandwidths. A scalar uses the same bandwidth for all dimensions, groups, and evaluation points. Otherwise, supply four values for each evaluation point: {it:h01 h02 h11 h12}, corresponding to the first and second running-variable coordinates in the control group and the first and second coordinates in the treated group. If {cmd:h()} is omitted, {cmd:rd2d} computes automatic bandwidths using {help rdbw2d:rdbw2d}.{p_end}

{p 4 8}{cmd:fuzzy(}{it:varname}{cmd:)} specifies the treatment receipt/status variable for fuzzy boundary discontinuity designs. The main table reports the fuzzy Wald ratio. The reduced-form outcome discontinuity is stored in {cmd:e(itt)}, and the first-stage treatment receipt/status discontinuity is stored in {cmd:e(fs)}. Side-specific fuzzy tables are stored in {cmd:e(itt_0)}, {cmd:e(itt_1)}, {cmd:e(fs_0)}, and {cmd:e(fs_1)}.{p_end}

{p 4 8}{cmd:deriv()} specifies derivative orders in the first and second coordinates. For example, {cmd:deriv(1 2)} estimates the mixed derivative with one derivative in the first coordinate and two in the second coordinate. Default is {cmd:deriv(0 0)}.{p_end}

{p 4 8}{cmd:tangvec()} supplies directional derivative vectors, one pair per boundary point. If supplied, {cmd:tangvec()} overrides {cmd:deriv()} and requires {cmd:p(1)} or higher.{p_end}

{dlgtab:Local Polynomial Regression}

{p 4 8}{cmd:p()} and {cmd:q()} set the local polynomial orders for point estimation and bias-corrected inference. Defaults are {cmd:p(1)} and {cmd:q(p+1)}.{p_end}

{p 4 8}{cmd:kernel()} specifies the kernel function. Available choices are:{p_end}
{p 8 12}{cmd:triangular} or {cmd:tri}: triangular kernel, the default.{p_end}
{p 8 12}{cmd:epanechnikov} or {cmd:epa}: Epanechnikov kernel.{p_end}
{p 8 12}{cmd:uniform} or {cmd:uni}: uniform kernel.{p_end}
{p 8 12}{cmd:gaussian} or {cmd:gau}: Gaussian kernel.{p_end}

{p 4 8}{cmd:kerneltype()} specifies the kernel structure. Available choices are:{p_end}
{p 8 12}{cmd:prod}: product kernel, the default.{p_end}
{p 8 12}{cmd:rad}: radial kernel.{p_end}

{dlgtab:Variance-Covariance Estimation}

{p 4 8}{cmd:vce()} specifies the variance-covariance estimator. Available choices are:{p_end}
{p 8 12}{cmd:hc0}: heteroskedasticity-robust plug-in residual variance estimator without small-sample adjustment.{p_end}
{p 8 12}{cmd:hc1}: heteroskedasticity-robust plug-in residual variance estimator with HC1 degrees-of-freedom adjustment, the default.{p_end}
{p 8 12}{cmd:hc2}: heteroskedasticity-robust plug-in residual variance estimator with HC2 leverage adjustment.{p_end}
{p 8 12}{cmd:hc3}: heteroskedasticity-robust plug-in residual variance estimator with HC3 leverage adjustment.{p_end}

{p 4 8}{cmd:cluster(}{it:clustvar}{cmd:)} specifies a cluster ID variable. When supplied, cluster-level influence sums with degrees-of-freedom weights are used.{p_end}

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
{p 8 12}{cmd:msetwo}: two different MSE-optimal bandwidth selectors, one for the control group and one for the treated group, at each evaluation point.{p_end}
{p 8 12}{cmd:certwo}: CER-optimal counterpart of {cmd:msetwo}.{p_end}
{p 8 12}{cmd:imsetwo}: two IMSE-optimal bandwidth selectors, one for each group, based on all evaluation points.{p_end}
{p 8 12}{cmd:icertwo}: CER-optimal counterpart of {cmd:imsetwo}.{p_end}
{p 8 12}{cmd:user provided}: used when {cmd:h()} is supplied.{p_end}

{p 4 8}{cmd:bwparam()} specifies the target parameter for fuzzy automatic bandwidth selection. Available choices are:{p_end}
{p 8 12}{cmd:main}: selects bandwidths for the linearized fuzzy Wald ratio, the default.{p_end}
{p 8 12}{cmd:itt}: selects bandwidths for the reduced-form outcome discontinuity. This option is ignored in sharp designs and when {cmd:h()} is supplied.{p_end}

{p 4 8}{cmd:method()} specifies the method for bandwidth calculations. Available choices are:{p_end}
{p 8 12}{cmd:dpi}: data-driven plug-in MSE-optimal selector, the default.{p_end}
{p 8 12}{cmd:rot}: rule-of-thumb selector.{p_end}

{p 4 8}{cmd:bwcheck()} enlarges preliminary bandwidths, when needed, so that at least the requested number of unique observations are included in each side-specific window. Default is {cmd:50 + p + 1}.{p_end}

{p 4 8}{cmd:masspoints()} specifies how repeated running-variable values are handled. Available choices are:{p_end}
{p 8 12}{cmd:check}: detects mass points and reports unique-observation counts, the default.{p_end}
{p 8 12}{cmd:adjust}: adjusts preliminary bandwidths to ensure a minimum number of unique observations within each side of the boundary.{p_end}
{p 8 12}{cmd:off}: ignores mass points.{p_end}

{p 4 8}{cmd:scaleregul()} specifies the scaling factor for the regularization term in bandwidth selection. Default is {cmd:scaleregul(3)} for {cmd:rd2d}.{p_end}

{p 4 8}{cmd:scalebiascrct()} specifies the scaling factor used for bias correction based on higher-order expansions. Default is {cmd:scalebiascrct(1)}.{p_end}

{p 4 8}{cmd:stdvars()} controls whether the running variables are standardized before automatic bandwidth selection. Use {cmd:stdvars(on)} or {cmd:stdvars(off)}. Default is {cmd:on}; this option affects automatic bandwidth selection only.{p_end}

{marker examples}{...}
{title:Examples}

{p 4 8}Simulated sharp design{p_end}
{p 8 8}{cmd:. clear}{p_end}
{p 8 8}{cmd:. set type double}{p_end}
{p 8 8}{cmd:. set obs 800}{p_end}
{p 8 8}{cmd:. set seed 123}{p_end}
{p 8 8}{cmd:. generate double x1 = rnormal()}{p_end}
{p 8 8}{cmd:. generate double x2 = rnormal()}{p_end}
{p 8 8}{cmd:. generate double d = x1 >= 0}{p_end}
{p 8 8}{cmd:. generate double y = 3 + 2*x1 + 1.5*x2 + d + rnormal()}{p_end}
{p 8 8}{cmd:. rd2d y x1 x2 d, b(0 0 0 1) h(.9) masspoints(off)}{p_end}

{p 4 8}Fuzzy design{p_end}
{p 8 8}{cmd:. generate double takeup = runiform() < cond(d, .8, .2)}{p_end}
{p 8 8}{cmd:. generate double yf = 3 + 2*x1 + 1.5*x2 + 1.5*takeup + rnormal()}{p_end}
{p 8 8}{cmd:. rd2d yf x1 x2 d, b(0 0 0 1) h(.9) fuzzy(takeup)}{p_end}

{p 4 8}Covariate adjustment and legacy fitting convention{p_end}
{p 8 8}{cmd:. generate double z = x1 + x2 + rnormal()}{p_end}
{p 8 8}{cmd:. rd2d y x1 x2 d, b(0 0 0 1) h(.9) covseff(z) fitmethod(joint)}{p_end}
{p 8 8}{cmd:. rd2d y x1 x2 d, b(0 0 0 1) h(.9) fitmethod(separate)}{p_end}

{marker stored_results}{...}
{title:Stored results}

{p 4 8}{cmd:rd2d} stores the following in {cmd:e()}:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(main)}}main treatment-effect table{p_end}
{synopt:{cmd:e(bw)}}bandwidth and effective sample-size table{p_end}
{synopt:{cmd:e(V_main)}}covariance matrix for the main bias-corrected estimates{p_end}
{synopt:{cmd:e(main_0)}}control-side table, sharp designs{p_end}
{synopt:{cmd:e(main_1)}}treated-side table, sharp designs{p_end}
{synopt:{cmd:e(itt)}}reduced-form table, fuzzy designs{p_end}
{synopt:{cmd:e(fs)}}first-stage table, fuzzy designs{p_end}
{synopt:{cmd:e(itt_0)}}control-side reduced-form table, fuzzy designs{p_end}
{synopt:{cmd:e(itt_1)}}treated-side reduced-form table, fuzzy designs{p_end}
{synopt:{cmd:e(fs_0)}}control-side first-stage table, fuzzy designs{p_end}
{synopt:{cmd:e(fs_1)}}treated-side first-stage table, fuzzy designs{p_end}
{synopt:{cmd:e(V_itt)}}covariance matrix for {cmd:e(itt)}{p_end}
{synopt:{cmd:e(V_fs)}}covariance matrix for {cmd:e(fs)}{p_end}
{synopt:{cmd:e(V_itt_0)}, {cmd:e(V_itt_1)}}covariance matrices for side-specific reduced-form tables{p_end}
{synopt:{cmd:e(V_fs_0)}, {cmd:e(V_fs_1)}}covariance matrices for side-specific first-stage tables{p_end}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(p)}}polynomial order for point estimation{p_end}
{synopt:{cmd:e(q)}}polynomial order for bias-corrected inference{p_end}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(fitmethod)}}fitting convention{p_end}
{synopt:{cmd:e(covseff)}}covariates used for efficiency adjustment, if supplied{p_end}

{marker references}{...}
{title:References}

{p 4 8}Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026). Estimation and Inference in Boundary Discontinuity Designs: Location-Based Methods.{p_end}

{p 4 8}Cattaneo, M. D., Titiunik, R., and Yu, R. R. (2026). Boundary Discontinuity Designs: Theory and Practice.{p_end}

{marker authors}{...}
{title:Authors}

{p 4 8}Matias D. Cattaneo, Princeton University. {browse "mailto:matias.d.cattaneo@gmail.com":matias.d.cattaneo@gmail.com}{p_end}
{p 4 8}Rocio Titiunik, Princeton University. {browse "mailto:rocio.titiunik@gmail.com":rocio.titiunik@gmail.com}{p_end}
{p 4 8}Ruiqi Rae Yu, Princeton University. {browse "mailto:raeyuuuu@gmail.com":raeyuuuu@gmail.com}{p_end}
