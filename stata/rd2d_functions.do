********************************************************************************
* RD2D STATA PACKAGE -- Mata functions
* Authors: Matias D. Cattaneo, Rocio Titiunik, Ruiqi Rae Yu
********************************************************************************
*!version 1.0.0  2026-05-26

version 16.0

mata:

struct rd2d_lfit {
	real rowvector estimate
	real rowvector se
	real matrix influence
	real scalar n_eff
}

struct rd2d_orderfit {
	real matrix mu0
	real matrix mu1
	real matrix se0
	real matrix se1
	real colvector N0
	real colvector N1
	real matrix infl0_1
	real matrix infl1_1
	real matrix infl0_2
	real matrix infl1_2
}

struct rd2d_bwfit {
	real colvector beta
	real matrix covconst
	real scalar eN
}

struct rd2d_bwconst {
	real scalar B
	real scalar V
	real scalar Reg1
	real scalar Reg2
}

struct rd2d_polyfit {
	real colvector coef
	real colvector se
}

struct rd2d_covcomp {
	real matrix zwz
	real matrix zwy
}

real scalar rd2d_mlib_loaded()
{
	return(1)
}

real matrix rd2d_cluster_sums(real matrix values, real colvector cluster, real colvector clusters)
{
	real matrix out, sorted_values, sums, info
	real colvector ord, sorted_cluster, uniq
	real scalar g, c

	out = J(rows(clusters), cols(values), 0)
	if (rows(cluster) == 0) return(out)

	ord = order(cluster, 1)
	sorted_cluster = cluster[ord]
	sorted_values = values[ord,.]
	info = panelsetup(sorted_cluster, 1)
	sums = panelsum(sorted_values, info)
	uniq = sorted_cluster[info[,1]]
	c = 1
	for (g=1; g<=rows(uniq); g++) {
		while (c <= rows(clusters) & clusters[c] < uniq[g]) c++
		if (c <= rows(clusters) & clusters[c] == uniq[g]) out[c,.] = sums[g,.]
	}
	return(out)
}

real matrix rd2d_parse_numlist(string scalar s, real scalar ncol)
{
	string rowvector toks
	real rowvector vals
	real matrix out
	real scalar i, j, n

	toks = tokens(strtrim(s))
	if (cols(toks) == 0) return(J(0, ncol, .))
	vals = strtoreal(toks)
	if (any(vals :>= .)) {
		errprintf("numeric option contains a nonnumeric value\n")
		_error(198)
	}
	if (mod(cols(vals), ncol) != 0) {
		errprintf("numeric option must contain a multiple of %g values\n", ncol)
		_error(198)
	}
	n = cols(vals) / ncol
	out = J(n, ncol, .)
	for (i=1; i<=n; i++) {
		for (j=1; j<=ncol; j++) out[i,j] = vals[(i-1)*ncol + j]
	}
	return(out)
}

real colvector rd2d_parse_kink_position(string scalar s, real scalar neval,
	real matrix b)
{
	string rowvector toks
	real rowvector vals
	real colvector out
	real scalar i, ind

	out = J(neval, 1, 0)
	s = strtrim(s)
	if (s == "") return(out)
	toks = tokens(s)
	vals = strtoreal(toks)
	if (any(vals :>= .)) {
		errprintf("kinkposition() contains a nonnumeric value\n")
		_error(198)
	}
	if (cols(vals) == neval & all((vals :== 0) :| (vals :== 1))) {
		out = vals'
	}
	else {
		for (i=1; i<=cols(vals); i++) {
			if (vals[i] != floor(vals[i]) | vals[i] < 1 | vals[i] > neval) {
				errprintf("kinkposition() must contain 1-based integer indices or a 0/1 indicator for each distance variable\n")
				_error(198)
			}
			ind = vals[i]
			out[ind] = 1
		}
	}
	if (sum(out) > 0 & any(b :>= .)) {
		errprintf("kinkposition() requires b() so distances between boundary points can be computed\n")
		_error(198)
	}
	return(out)
}

real colvector rd2d_distance_to_known_kink(real matrix b, real colvector kinkpos)
{
	real colvector idx, out, d
	real scalar j

	out = J(rows(b), 1, .)
	if (sum(kinkpos) == 0) return(out)
	idx = selectindex(kinkpos :== 1)
	for (j=1; j<=rows(b); j++) {
		d = sqrt((b[j,1] :- b[idx,1]):^2 + (b[j,2] :- b[idx,2]):^2)
		out[j] = min(d)
	}
	return(out)
}

real scalar rd2d_fact(real scalar n)
{
	real scalar out, i
	out = 1
	if (n <= 1) return(out)
	for (i=2; i<=n; i++) out = out * i
	return(out)
}

real scalar rd2d_cer_factor(real scalar n, real scalar p)
{
	n = max((n, 1))
	return(n^(1/(2*p+4) - 1/(p+4)))
}

string scalar rd2d_bwbase(string scalar bwselect)
{
	bwselect = strlower(strtrim(bwselect))
	if (bwselect == "cerrd") return("mserd")
	if (bwselect == "certwo") return("msetwo")
	if (bwselect == "icerrd") return("imserd")
	if (bwselect == "icertwo") return("imsetwo")
	return(bwselect)
}

real scalar rd2d_is_cer(string scalar bwselect)
{
	bwselect = strlower(strtrim(bwselect))
	return(anyof(("cerrd","certwo","icerrd","icertwo"), bwselect))
}

real scalar rd2d_is_common(string scalar bwselect)
{
	string scalar base
	base = rd2d_bwbase(bwselect)
	return(anyof(("mserd","imserd"), base))
}

real matrix rd2d_kernel(real matrix u, string scalar kernel)
{
	real matrix au, w
	kernel = strlower(strtrim(kernel))
	au = abs(u)
	if (kernel == "tri" | kernel == "triangular" | kernel == "") {
		w = (1 :- au) :* (au :<= 1)
	}
	else if (kernel == "epa" | kernel == "epanechnikov") {
		w = .75 :* (1 :- u:^2) :* (au :<= 1)
	}
	else if (kernel == "uni" | kernel == "uniform") {
		w = .5 :* (au :<= 1)
	}
	else if (kernel == "gau" | kernel == "gaussian") {
		w = exp(-.5 :* u:^2) :/ sqrt(2*pi())
	}
	else {
		errprintf("kernel() incorrectly specified\n")
		_error(198)
	}
	return(w)
}

real colvector rd2d_kernel_nonnegative(real colvector u, string scalar kernel)
{
	real colvector w
	kernel = strlower(strtrim(kernel))
	if (kernel == "tri" | kernel == "triangular" | kernel == "") {
		w = (1 :- u) :* (u :<= 1)
	}
	else if (kernel == "epa" | kernel == "epanechnikov") {
		w = .75 :* (1 :- u:^2) :* (u :<= 1)
	}
	else if (kernel == "uni" | kernel == "uniform") {
		w = .5 :* (u :<= 1)
	}
	else if (kernel == "gau" | kernel == "gaussian") {
		w = exp(-.5 :* u:^2) :/ sqrt(2*pi())
	}
	else {
		w = rd2d_kernel(u, kernel)
	}
	return(w)
}

real matrix rd2d_basis1(real colvector x, real scalar p)
{
	real matrix R
	real scalar j
	R = J(rows(x), p+1, 1)
	for (j=2; j<=p+1; j++) R[,j] = R[,j-1] :* x
	return(R)
}

real matrix rd2d_basis2(real matrix x, real scalar p)
{
	real matrix R
	real scalar n, k, deg, ypow, xpow, col
	n = rows(x)
	k = (p+1)*(p+2)/2
	R = J(n, k, 1)
	col = 1
	for (deg=0; deg<=p; deg++) {
		for (ypow=0; ypow<=deg; ypow++) {
			xpow = deg - ypow
			R[,col] = (x[,1]:^xpow) :* (x[,2]:^ypow)
			col++
		}
	}
	return(R)
}

real rowvector rd2d_target2(real scalar p, real rowvector deriv, real matrix tangvec, real scalar j)
{
	real rowvector target
	real scalar deg, ypow, xpow, col
	target = J(1, (p+1)*(p+2)/2, 0)
	if (rows(tangvec) > 0) {
		if (p < 1) {
			errprintf("tangvec() requires p >= 1\n")
			_error(198)
		}
		target[2] = tangvec[j,1]
		target[3] = tangvec[j,2]
		return(target)
	}
	col = 1
	for (deg=0; deg<=p; deg++) {
		for (ypow=0; ypow<=deg; ypow++) {
			xpow = deg - ypow
			if (xpow == deriv[1] & ypow == deriv[2]) {
				target[col] = rd2d_fact(deriv[1]) * rd2d_fact(deriv[2])
				return(target)
			}
			col++
		}
	}
	return(target)
}

real rowvector rd2d_target1(real scalar p)
{
	real rowvector target
	target = J(1, p+1, 0)
	target[1] = 1
	return(target)
}

real colvector rd2d_location_weights(real matrix centered, real rowvector h, string scalar kernel, string scalar kerneltype)
{
	real colvector w, u
	real scalar hprod
	hprod = max((h[1]*h[2], epsilon(1)))
	kerneltype = strlower(strtrim(kerneltype))
	if (kerneltype == "rad") {
		u = sqrt((centered[,1]:/h[1]):^2 + (centered[,2]:/h[2]):^2)
		w = rd2d_kernel(u, kernel) :/ hprod
	}
	else {
		w = (rd2d_kernel(centered[,1]:/h[1], kernel) :*
		     rd2d_kernel(centered[,2]:/h[2], kernel)) :/ hprod
	}
	return(w)
}

real colvector rd2d_distance_weights(real colvector sdist, real scalar h, real scalar treated, string scalar kernel)
{
	real colvector u, w
	u = abs(sdist)
	w = rd2d_kernel(u:/h, kernel) :/ max((h^2, epsilon(1)))
	if (treated) w = w :* (sdist :>= 0)
	else         w = w :* (sdist :<  0)
	return(w)
}

struct rd2d_lfit scalar rd2d_local_fit(real matrix design, real colvector weights,
	real matrix outcomes, real rowvector target, string scalar vce, real colvector cluster)
{
	return(rd2d_local_fit_scaled(design, weights, outcomes, target, vce, cluster, ., J(0,1,.)))
}

struct rd2d_lfit scalar rd2d_local_fit_scaled(real matrix design, real colvector weights,
	real matrix outcomes, real rowvector target, string scalar vce, real colvector cluster,
	real scalar scale_override, real colvector cluster_groups)
{
	struct rd2d_lfit scalar out
	real colvector idx, w, leverage, adj, score, sqrtadj, ckeep, groups_eff, all_groups
	real matrix X, Y, gram, invG, beta, fitted, resid, Xw, infl_kept, infl_full, cov
	real rowvector row
	real scalar n, k, m, i, scale, has_cluster

	n = rows(design)
	k = cols(design)
	m = cols(outcomes)
	has_cluster = (rows(cluster) == n)
	out.estimate = J(1, m, .)
	out.se = J(1, m, .)
	out.n_eff = 0
	all_groups = J(0,1,.)
	if (has_cluster) {
		all_groups = rows(cluster_groups) > 0 ? cluster_groups : uniqrows(sort(cluster,1))
	}
	out.influence = J(has_cluster ? rows(all_groups) : n, m, 0)

	idx = selectindex(weights :> 0)
	if (rows(idx) == 0) return(out)
	X = design[idx,.]
	Y = outcomes[idx,.]
	w = weights[idx]
	n = rows(X)
	k = cols(X)
	out.n_eff = n

	Xw = X :* (w * J(1,k,1))
	gram = X' * Xw
	invG = pinv(gram)
	beta = invG * (X' * (Y :* (w * J(1,m,1))))
	fitted = X * beta
	resid = Y - fitted

	adj = J(n,1,1)
	vce = strlower(strtrim(vce))
	if (vce == "hc2" | vce == "hc3") {
		leverage = rowsum((X * invG) :* X) :* w
		if (vce == "hc2") adj = 1 :/ rowmax((1 :- leverage, J(n,1,1e-8)))
		if (vce == "hc3") adj = (1 :/ rowmax((1 :- leverage, J(n,1,1e-8)))):^2
	}

	row = target * invG
	score = Xw * row'
	sqrtadj = sqrt(adj)
	infl_kept = (score * J(1,m,1)) :* resid :* (sqrtadj * J(1,m,1))

	scale = 1
	if (scale_override < .) scale = scale_override
	if (has_cluster) {
		ckeep = cluster[idx]
		if (scale_override >= . & n > k) groups_eff = uniqrows(sort(ckeep,1))
		infl_full = rd2d_cluster_sums(infl_kept, ckeep, all_groups)
		if (scale_override >= . & n > k) {
			scale = 1
			if (vce == "hc1") scale = scale * n/(n-k)
			if (rows(groups_eff) > 1) {
				scale = scale * ((n-1)/(n-k)) * (rows(groups_eff)/(rows(groups_eff)-1))
			}
		}
		cov = scale * infl_full' * infl_full
		out.influence = sqrt(scale) * infl_full
	}
	else {
		if (scale_override >= . & vce == "hc1" & n > k) scale = n/(n-k)
		cov = scale * infl_kept' * infl_kept
		infl_full = J(rows(design), m, 0)
		infl_full[idx,.] = sqrt(scale) * infl_kept
		out.influence = infl_full
	}
	out.estimate = target * beta
	out.se = sqrt(diagonal(cov))'
	return(out)
}

real scalar rd2d_joint_vce_scale(string scalar vce, real scalar eN,
	real scalar kdf, real colvector cluster, real scalar clustered)
{
	real scalar scale, g
	scale = 1
	vce = strlower(strtrim(vce))
	if (eN <= kdf) return(scale)
	if (vce == "hc1") scale = scale * eN / (eN - kdf)
	if (clustered & rows(cluster) > 0) {
		g = rows(uniqrows(sort(cluster,1)))
		if (g > 1) scale = scale * ((eN - 1)/(eN - kdf)) * (g/(g - 1))
	}
	return(scale)
}

struct rd2d_covcomp scalar rd2d_projection_components(real matrix design,
	real colvector weights, real matrix outcomes, real matrix covs)
{
	struct rd2d_covcomp scalar out
	real colvector idx, w
	real matrix X, Y, Z, invG, Uy, Uz
	real scalar zcols

	zcols = cols(covs)
	out.zwz = J(zcols, zcols, 0)
	out.zwy = J(zcols, cols(outcomes), 0)
	idx = selectindex(weights :> 0)
	if (rows(idx) == 0 | zcols == 0) return(out)
	X = design[idx,.]
	Y = outcomes[idx,.]
	Z = covs[idx,.]
	w = weights[idx]
	invG = pinv(X' * (X :* (w * J(1, cols(X), 1))))
	Uy = X' * (Y :* (w * J(1, cols(Y), 1)))
	Uz = X' * (Z :* (w * J(1, cols(Z), 1)))
	out.zwz = Z' * (Z :* (w * J(1, cols(Z), 1))) - Uz' * invG * Uz
	out.zwy = Z' * (Y :* (w * J(1, cols(Y), 1))) - Uz' * invG * Uy
	return(out)
}

real matrix rd2d_covariate_gamma(real matrix design, real colvector w0,
	real colvector w1, real matrix outcomes, real matrix covs)
{
	struct rd2d_covcomp scalar c0, c1
	if (cols(covs) == 0) return(J(0, cols(outcomes), .))
	c0 = rd2d_projection_components(design, w0, outcomes, covs)
	c1 = rd2d_projection_components(design, w1, outcomes, covs)
	return(pinv(c0.zwz + c1.zwz) * (c0.zwy + c1.zwy))
}

real scalar rd2d_covariate_rank(real matrix design, real colvector w0,
	real colvector w1, real matrix covs)
{
	struct rd2d_covcomp scalar c0, c1
	if (cols(covs) == 0) return(0)
	c0 = rd2d_projection_components(design, w0, J(rows(design),1,0), covs)
	c1 = rd2d_projection_components(design, w1, J(rows(design),1,0), covs)
	return(rank(c0.zwz + c1.zwz))
}

real matrix rd2d_apply_covariates(real matrix outcomes, real matrix covs,
	real matrix gamma)
{
	if (cols(covs) == 0 | rows(gamma) == 0) return(outcomes)
	return(outcomes - covs * gamma)
}

real rowvector rd2d_bw_vce_scales(string scalar fitmethod, string scalar vce,
	real scalar covadj, real scalar e0, real scalar e1, real scalar kdf,
	real scalar n_cov)
{
	real rowvector out
	real scalar scale

	out = (., .)
	fitmethod = strlower(strtrim(fitmethod))
	vce = strlower(strtrim(vce))
	if (vce != "hc1" | (fitmethod != "joint" & !covadj)) return(out)
	if (fitmethod == "joint") {
		scale = rd2d_joint_vce_scale(vce, e0 + e1, 2*kdf + n_cov, J(0,1,.), 0)
		return((scale, scale))
	}
	return((rd2d_joint_vce_scale(vce, e0, kdf + n_cov, J(0,1,.), 0),
		rd2d_joint_vce_scale(vce, e1, kdf + n_cov, J(0,1,.), 0)))
}

real scalar rd2d_basis_count2(real scalar p)
{
	return((p+1)*(p+2)/2)
}

real matrix rd2d_H2(real scalar h, real scalar p)
{
	real colvector dvec
	real scalar deg, ypow, xpow, col
	dvec = J(rd2d_basis_count2(p), 1, 1)
	col = 1
	for (deg=0; deg<=p; deg++) {
		for (ypow=0; ypow<=deg; ypow++) {
			xpow = deg - ypow
			dvec[col] = h^(xpow + ypow)
			col++
		}
	}
	return(diag(dvec))
}

real matrix rd2d_invH2(real scalar h, real scalar p)
{
	real matrix H
	H = rd2d_H2(h, p)
	return(diag(1 :/ diagonal(H)))
}

real matrix rd2d_vce_const2(real matrix wR, real colvector resd, real scalar h,
	real colvector cluster)
{
	real matrix scores, summed, part
	real colvector groups
	real scalar g, n, k, scale

	scores = wR :* (resd * J(1, cols(wR), 1))
	if (rows(cluster) == rows(resd)) {
		groups = uniqrows(sort(cluster, 1))
		summed = rd2d_cluster_sums(scores, cluster, groups)
		n = rows(cluster)
		k = cols(scores)
		scale = (rows(groups) > 1 & n > k) ? ((n - 1) / (n - k)) * (rows(groups) / (rows(groups) - 1)) : 1
		return(scale * summed' * summed * h^2)
	}
	return(scores' * scores * h^2)
}

struct rd2d_bwfit scalar rd2d_lm_exact2(real matrix x, real colvector y,
	real colvector distance, real scalar h, real scalar p, string scalar vce,
	string scalar kernel, string scalar kerneltype, real colvector cluster,
	real scalar varr)
{
	struct rd2d_bwfit scalar out
	real colvector w, idx, ew, eY, sqrtw, resd, hii, ecluster
	real matrix eu, eR, sqrtw_R, invG, beta, H, invH, cov
	real scalar n, k

	out.beta = J(rd2d_basis_count2(p), 1, .)
	out.covconst = J(rd2d_basis_count2(p), rd2d_basis_count2(p), .)
	out.eN = 0
	w = rd2d_location_weights(x, (h, h), kernel, kerneltype)
	idx = selectindex(w :> 0)
	if (rows(idx) == 0) return(out)

	ew = w[idx]
	eY = y[idx]
	eu = x[idx,.] :/ h
	eR = rd2d_basis2(eu, p)
	sqrtw = sqrt(ew)
	sqrtw_R = eR :* (sqrtw * J(1, cols(eR), 1))
	invG = pinv(sqrtw_R' * sqrtw_R)
	H = rd2d_H2(h, p)
	invH = rd2d_invH2(h, p)
	beta = invH * invG * (sqrtw_R' * (sqrtw :* eY))
	out.beta = beta
	out.eN = rows(idx)

	if (varr) {
		resd = eY - eR * (H * beta)
		n = rows(eR)
		k = cols(eR)
		vce = strlower(strtrim(vce))
		if (vce == "hc1" & n > k) {
			resd = resd :* sqrt(n / (n - k))
		}
		else if (vce == "hc2" | vce == "hc3") {
			hii = rowsum((sqrtw_R * invG) :* sqrtw_R)
			if (vce == "hc2") resd = resd :* sqrt(1 :/ rowmax((1 :- hii, J(n,1,1e-12))))
			else              resd = resd :* (1 :/ rowmax((1 :- hii, J(n,1,1e-12))))
		}
		ecluster = rows(cluster) == rows(y) ? cluster[idx] : J(0,1,.)
		cov = rd2d_vce_const2(eR :* (ew * J(1, cols(eR), 1)), resd, h, ecluster)
		out.covconst = invG' * cov * invG
	}
	return(out)
}

real matrix rd2d_lm_multi_exact2(real matrix x, real matrix outcomes,
	real colvector distance, real scalar h, real scalar p, string scalar kernel,
	string scalar kerneltype)
{
	real colvector w, idx, ew, sqrtw
	real matrix eu, eR, sqrtw_R, invG, H, invH, eY

	w = rd2d_location_weights(x, (h, h), kernel, kerneltype)
	idx = selectindex(w :> 0)
	if (rows(idx) == 0) return(J(rd2d_basis_count2(p), cols(outcomes), .))
	ew = w[idx]
	eu = x[idx,.] :/ h
	eR = rd2d_basis2(eu, p)
	sqrtw = sqrt(ew)
	sqrtw_R = eR :* (sqrtw * J(1, cols(eR), 1))
	invG = pinv(sqrtw_R' * sqrtw_R)
	H = rd2d_H2(h, p)
	invH = rd2d_invH2(h, p)
	eY = outcomes[idx,.]
	return(invH * invG * (sqrtw_R' * (eY :* (sqrtw * J(1, cols(eY), 1)))))
}

real rowvector rd2d_get_coeff_exact2(real matrix x, real colvector distance,
	real rowvector vec, real scalar p, real scalar h, string scalar kernel,
	string scalar kerneltype)
{
	real colvector w, idx, ew, sqrtw
	real matrix eu, eRaug, eR, eS, sqrtw_R, sqrtw_S, invG
	real scalar pcount, p1count
	real rowvector tail

	w = rd2d_location_weights(x, (h, h), kernel, kerneltype)
	idx = selectindex(w :> 0)
	pcount = rd2d_basis_count2(p)
	p1count = rd2d_basis_count2(p + 1)
	if (rows(idx) == 0) return(J(1, p1count, 0))
	ew = w[idx]
	eu = x[idx,.] :/ h
	eRaug = rd2d_basis2(eu, p + 1)
	eR = eRaug[,1..pcount]
	eS = eRaug[,(pcount+1)..p1count]
	sqrtw = sqrt(ew)
	sqrtw_R = eR :* (sqrtw * J(1, cols(eR), 1))
	sqrtw_S = eS :* (sqrtw * J(1, cols(eS), 1))
	invG = pinv(sqrtw_R' * sqrtw_R)
	tail = vec * invG * sqrtw_R' * sqrtw_S
	return(J(1, pcount, 0), tail)
}

struct rd2d_bwconst scalar rd2d_bw_v2_exact2_adj(real matrix x, real colvector y,
	real colvector distance, real scalar p, real rowvector vec, real scalar dn,
	real scalar bn1, real scalar bn2, string scalar vce, string scalar kernel,
	string scalar kerneltype, real colvector cluster, real colvector yvar,
	real colvector ybias, real colvector yreg2, real scalar scale_var,
	real scalar scale_bias)
{
	struct rd2d_bwconst scalar out
	struct rd2d_bwfit scalar fitp1, fitp2
	real colvector w, idx, ew, eY, sqrtw, sqrtw_Y, resd, hii, ecluster
	real matrix eu, eRaug, eR, eS, eT, sqrtw_R, sqrtw_S, sqrtw_T, invG
	real matrix H, invH, beta, sigma
	real rowvector vecq, vect
	real scalar pcount, p1count, p2count, n, k
	string scalar vce_bias

	out.B = .
	out.V = .
	out.Reg1 = .
	out.Reg2 = .
	if (rows(yvar) != rows(y)) yvar = y
	if (rows(ybias) != rows(y)) ybias = y
	if (rows(yreg2) != rows(y)) yreg2 = y
	w = rd2d_location_weights(x, (dn, dn), kernel, kerneltype)
	idx = selectindex(w :> 0)
	if (rows(idx) == 0) return(out)
	ew = w[idx]
	eY = yvar[idx]
	eu = x[idx,.] :/ dn
	pcount = rd2d_basis_count2(p)
	p1count = rd2d_basis_count2(p + 1)

	if (bn2 >= .) {
		eRaug = rd2d_basis2(eu, p + 1)
		eR = eRaug[,1..pcount]
		eS = eRaug[,(pcount+1)..p1count]
		eT = J(rows(eR), 0, .)
	}
	else {
		p2count = rd2d_basis_count2(p + 2)
		eRaug = rd2d_basis2(eu, p + 2)
		eR = eRaug[,1..pcount]
		eS = eRaug[,(pcount+1)..p1count]
		eT = eRaug[,(p1count+1)..p2count]
	}

	sqrtw = sqrt(ew)
	sqrtw_R = eR :* (sqrtw * J(1, cols(eR), 1))
	sqrtw_S = eS :* (sqrtw * J(1, cols(eS), 1))
	sqrtw_Y = sqrtw :* eY
	invG = pinv(sqrtw_R' * sqrtw_R)
	vecq = J(1, pcount, 0), vec * invG * sqrtw_R' * sqrtw_S
	vect = J(1, 0, .)
	if (cols(eT) > 0) {
		sqrtw_T = eT :* (sqrtw * J(1, cols(eT), 1))
		vect = J(1, p1count, 0), vec * invG * sqrtw_R' * sqrtw_T
	}

	H = rd2d_H2(dn, p)
	invH = rd2d_invH2(dn, p)
	beta = invH * invG * (sqrtw_R' * sqrtw_Y)
	resd = eY - eR * (H * beta)
	n = rows(eR)
	k = cols(eR)
	vce = strlower(strtrim(vce))
	if (scale_var >= . & vce == "hc1" & n > k) {
		resd = resd :* sqrt(n / (n - k))
	}
	else if (vce == "hc2" | vce == "hc3") {
		hii = rowsum((sqrtw_R * invG) :* sqrtw_R)
		if (vce == "hc2") resd = resd :* sqrt(1 :/ rowmax((1 :- hii, J(n,1,1e-12))))
		else              resd = resd :* (1 :/ rowmax((1 :- hii, J(n,1,1e-12))))
	}
	ecluster = rows(cluster) == rows(y) ? cluster[idx] : J(0,1,.)
	sigma = rd2d_vce_const2(eR :* (ew * J(1, cols(eR), 1)), resd, dn, ecluster)
	out.V = (vec * invG' * sigma * invG * vec')[1,1]
	if (scale_var < .) out.V = out.V * scale_var

	vce_bias = (scale_bias < . & vce == "hc1") ? "hc0" : vce
	fitp1 = rd2d_lm_exact2(x, ybias, distance, bn1, p + 1, vce_bias, kernel, kerneltype, cluster, 1)
	out.B = (vecq * fitp1.beta)[1,1]
	out.Reg1 = (vecq * fitp1.covconst * vecq' / (bn1^(2 + 2*(p + 1))))[1,1]
	if (scale_bias < .) out.Reg1 = out.Reg1 * scale_bias
	if (bn2 < . & cols(vect) > 0) {
		fitp2 = rd2d_lm_exact2(x, yreg2, distance, bn2, p + 2, vce, kernel, kerneltype, cluster, 0)
		out.Reg2 = ((dn * vect) * fitp2.beta)[1,1]
	}
	return(out)
}

struct rd2d_bwconst scalar rd2d_bw_v2_exact2(real matrix x, real colvector y,
	real colvector distance, real scalar p, real rowvector vec, real scalar dn,
	real scalar bn1, real scalar bn2, string scalar vce, string scalar kernel,
	string scalar kerneltype, real colvector cluster)
{
	return(rd2d_bw_v2_exact2_adj(x, y, distance, p, vec, dn, bn1, bn2,
		vce, kernel, kerneltype, cluster, y, y, y, ., .))
}

real scalar rd2d_rot_location(real matrix x, string scalar kernel, real scalar M)
{
	real matrix cov, invcov
	real scalar detcov, sqrtdet, traceconst, mu2, l2, n_eff
	kernel = strlower(strtrim(kernel))
	if (kernel == "epa" | kernel == "epanechnikov") {
		mu2 = 1/6
		l2 = 4/(3*pi())
	}
	else if (kernel == "uni" | kernel == "uniform") {
		mu2 = 1/4
		l2 = 1/pi()
	}
	else if (kernel == "gau" | kernel == "gaussian") {
		mu2 = 1
		l2 = 1/(4*pi())
	}
	else {
		mu2 = 3/20
		l2 = 3/(2*pi())
	}
	cov = variance(x)
	invcov = pinv(cov)
	detcov = det(cov)
	if (detcov <= 0 | detcov >= .) detcov = max((epsilon(1), cov[1,1]*cov[2,2]))
	sqrtdet = sqrt(detcov)
	traceconst = (1/(2^4*pi()*sqrtdet)) * (2*trace(invcov*invcov) + trace(invcov)^2)
	n_eff = max((M, 1))
	return(((2*l2)/(n_eff*mu2*traceconst))^(1/6))
}

real scalar rd2d_rot_distance(real colvector u, string scalar kernel)
{
	real scalar mu2, l2, varhat, traceconst
	kernel = strlower(strtrim(kernel))
	if (kernel == "epa" | kernel == "epanechnikov") {
		mu2 = 1/6
		l2 = 4/(3*pi())
	}
	else if (kernel == "uni" | kernel == "uniform") {
		mu2 = 1/4
		l2 = 1/pi()
	}
	else if (kernel == "gau" | kernel == "gaussian") {
		mu2 = 1
		l2 = 1/(4*pi())
	}
	else {
		mu2 = 3/20
		l2 = 3/(2*pi())
	}
	varhat = .5 * mean(u:^2)
	if (varhat <= 0 | varhat >= .) varhat = variance(u)
	if (varhat <= 0 | varhat >= .) varhat = 1
	traceconst = 1/(2*pi()*varhat^3)
	return(((2*l2)/(rows(u)*mu2*traceconst))^(1/6))
}

real rowvector rd2d_kth_bounds(real colvector vals, real scalar bwcheck)
{
	real colvector s
	real scalar n, k
	s = sort(vals, 1)
	n = rows(s)
	if (n == 0) return((., .))
	k = min((max((bwcheck,1)), n))
	return((s[k], s[n]))
}

real scalar rd2d_quantile(real colvector vals, real scalar prob)
{
	real colvector s
	real scalar n, h, lo, hi

	s = sort(vals, 1)
	n = rows(s)
	if (n == 0) return(.)
	if (prob <= 0) return(s[1])
	if (prob >= 1) return(s[n])
	h = (n - 1) * prob + 1
	lo = floor(h)
	hi = ceil(h)
	if (lo == hi) return(s[lo])
	return(s[lo] + (h - lo) * (s[hi] - s[lo]))
}

real matrix rd2d_apply_location_bwcheck(real matrix H, real matrix x, real colvector d,
	real matrix b, real scalar bwcheck, string scalar kerneltype)
{
	real matrix out, centered
	real colvector dist, side0, side1, sd
	real rowvector b0, b1
	real scalar j

	out = H
	if (bwcheck <= 0) return(out)
	sd = sqrt(diagonal(variance(x)))
	sd = editmissing(sd, 1)
	sd = sd :+ (sd :<= 0)
	side0 = (d :== 0)
	side1 = (d :!= 0)
	for (j=1; j<=rows(b); j++) {
		centered = x :- (J(rows(x),1,1) * b[j,.])
		if (strlower(strtrim(kerneltype)) == "rad") {
			dist = sqrt(centered[,1]:^2 + centered[,2]:^2)
			b0 = rd2d_kth_bounds(select(dist, side0), bwcheck)
			b1 = rd2d_kth_bounds(select(dist, side1), bwcheck)
			out[j,1] = min((max((out[j,1], b0[1])), b0[2]))
			out[j,2] = min((max((out[j,2], b0[1])), b0[2]))
			out[j,3] = min((max((out[j,3], b1[1])), b1[2]))
			out[j,4] = min((max((out[j,4], b1[1])), b1[2]))
		}
		else {
			dist = rowmax((abs(centered[,1]:/sd[1]), abs(centered[,2]:/sd[2])))
			b0 = rd2d_kth_bounds(select(dist, side0), bwcheck)
			b1 = rd2d_kth_bounds(select(dist, side1), bwcheck)
			out[j,1] = min((max((out[j,1], b0[1]*sd[1])), b0[2]*sd[1]))
			out[j,2] = min((max((out[j,2], b0[1]*sd[2])), b0[2]*sd[2]))
			out[j,3] = min((max((out[j,3], b1[1]*sd[1])), b1[2]*sd[1]))
			out[j,4] = min((max((out[j,4], b1[1]*sd[2])), b1[2]*sd[2]))
		}
	}
	return(out)
}

real matrix rd2d_location_bw_simple(real colvector y, real matrix x, real colvector d,
	real matrix b, string scalar hstr, real scalar p, string scalar kernel,
	string scalar kerneltype, string scalar bwselect, real scalar bwcheck,
	real scalar stdvars)
{
	real matrix H, hraw, xw, bw
	real colvector sd
	real scalar hrot, j, neval, n0, n1, M, factor

	neval = rows(b)
	hraw = rd2d_parse_numlist(hstr, 1)
	if (rows(hraw) == 1 & cols(hraw) == 1) {
		H = J(neval, 4, hraw[1,1])
	}
	else if (rows(hraw) > 0) {
		if (rows(hraw)*cols(hraw) != neval*4) {
			errprintf("h() must be a scalar or contain 4 values per evaluation point\n")
			_error(198)
		}
		H = rd2d_parse_numlist(hstr, 4)
	}
	else {
		xw = x
		sd = J(2,1,1)
		if (stdvars) {
			sd = sqrt(diagonal(variance(x)))
			sd = editmissing(sd, 1)
			sd = sd :+ (sd :<= 0)
			xw = (x[,1]:/sd[1], x[,2]:/sd[2])
		}
		M = rows(uniqrows((xw,d)))
		hrot = rd2d_rot_location(xw, kernel, M)
		H = J(neval, 4, .)
		for (j=1; j<=neval; j++) {
			H[j,.] = (hrot*sd[1], hrot*sd[2], hrot*sd[1], hrot*sd[2])
		}
		if (rd2d_is_cer(bwselect)) {
			if (rd2d_is_common(bwselect)) {
				factor = rd2d_cer_factor(M, p)
				H = H :* factor
			}
			else {
				n0 = sum(d:==0)
				n1 = sum(d:!=0)
				H[,(1,2)] = H[,(1,2)] :* rd2d_cer_factor(n0, p)
				H[,(3,4)] = H[,(3,4)] :* rd2d_cer_factor(n1, p)
			}
		}
	}
	H = rd2d_apply_location_bwcheck(H, x, d, b, bwcheck, kerneltype)
	bw = b, H, J(neval,1,0), J(neval,1,0)
	for (j=1; j<=neval; j++) {
		bw[j,7] = sum(rd2d_location_weights(x :- (J(rows(x),1,1)*b[j,.]), H[j,(1,2)], kernel, kerneltype) :* (d:==0) :> 0)
		bw[j,8] = sum(rd2d_location_weights(x :- (J(rows(x),1,1)*b[j,.]), H[j,(3,4)], kernel, kerneltype) :* (d:!=0) :> 0)
	}
	return(bw)
}

real matrix rd2d_location_bw_full(real colvector y, real matrix x, real colvector d,
	real matrix b, string scalar hstr, real colvector fuzzy, real colvector cluster,
	real scalar p, real rowvector deriv, real matrix tangvec, string scalar kernel,
	string scalar kerneltype, string scalar bwselect, string scalar bwparam,
	string scalar method, string scalar vce, real scalar bwcheck,
	string scalar masspoints, real scalar scaleregul, real scalar scalebiascrct,
	real scalar stdvars, string scalar fitmethod, real matrix covs)
{
	real matrix hraw, xw, bwork, bw, results, H, beta0, beta1, outcomes
	real matrix design, gamma, adjusted, w0full, w1full
	real colvector sd, centered_dist, dist0, dist1, ybw, side0, side1
	real colvector w0dn, w1dn, w0thr, w1thr, w0bn, w1bn
	real colvector ybase, fuzzybase, ybnv, ybnb, yhnv, yhnb, yhnt
	real colvector ybw0, ybw1, ybnv0, ybnv1, ybnb0, ybnb1
	real colvector yhnv0, yhnv1, yhnb0, yhnb1, yhnt0, yhnt1
	real colvector cluster0, cluster1
	real matrix centered, x0, x1
	real rowvector target, vecq0, vecq1, b0, b1, scales_v, scales_b
	struct rd2d_bwconst scalar bnconst0, bnconst1, hnconst0, hnconst1
	real scalar neval, n, j, dn, dn0, dn1, bwmin0, bwmin1, bwmax0, bwmax1
	real scalar M, M0, M1, derivsum, derivdenom, thr0, thr1, bn0, bn1
	real scalar hn0, hn1, hn, tauitt, taufs, graditt, gradfs, i
	real scalar VV, BB, BB0, BB1, factor0, factor1, covadj, n_cov, clustered
	real scalar eN0, eN1, eNdn0, eNdn1, eNthr0, eNthr1, eNbn0, eNbn1

	neval = rows(b)
	hraw = rd2d_parse_numlist(hstr, 1)
	if (rows(hraw) > 0) {
		return(rd2d_location_bw_simple(y, x, d, b, hstr, p, kernel,
			kerneltype, bwselect, bwcheck, stdvars))
	}

	n = rows(x)
	xw = x
	bwork = b
	sd = J(2,1,1)
	if (stdvars) {
		sd = sqrt(diagonal(variance(x)))
		sd = editmissing(sd, 1)
		sd = sd :+ (sd :<= 0)
		xw = (x[,1]:/sd[1], x[,2]:/sd[2])
		bwork = (b[,1]:/sd[1], b[,2]:/sd[2])
	}

	M0 = sum(d :== 0)
	M1 = sum(d :!= 0)
	M = M0 + M1
	if (strlower(strtrim(masspoints)) == "check" | strlower(strtrim(masspoints)) == "adjust") {
		M0 = rows(uniqrows(select(xw,d:==0)))
		M1 = rows(uniqrows(select(xw,d:!=0)))
		M = M0 + M1
	}

	dn = rd2d_rot_location(xw, kernel, M)
	target = rd2d_target2(p, deriv, tangvec, 1)
	derivsum = sum(deriv)
	if (rows(tangvec) > 0) derivsum = 1
	derivdenom = 2*p + 2 - 2*derivsum
	results = J(neval, 16, .)
	bw = J(neval, 8, .)
	side0 = d :== 0
	side1 = d :!= 0
	clustered = (rows(cluster) == n)
	cluster0 = clustered ? select(cluster, side0) : J(0,1,.)
	cluster1 = clustered ? select(cluster, side1) : J(0,1,.)
	method = strlower(strtrim(method))
	bwselect = strlower(strtrim(bwselect))
	bwparam = strlower(strtrim(bwparam))
	fitmethod = (strlower(strtrim(fitmethod)) == "separate" ? "separate" : "joint")
	covadj = (cols(covs) > 0)

	for (j=1; j<=neval; j++) {
		target = rd2d_target2(p, deriv, tangvec, j)
		centered = xw :- (J(n,1,1) * bwork[j,.])
		centered_dist = sqrt(centered[,1]:^2 + centered[,2]:^2)
		dist0 = select(centered_dist, side0)
		dist1 = select(centered_dist, side1)
		x0 = select(centered, side0)
		x1 = select(centered, side1)
		dn0 = dn
		dn1 = dn
		bwmin0 = 0
		bwmin1 = 0
		bwmax0 = .
		bwmax1 = .
		if (bwcheck > 0) {
			b0 = rd2d_kth_bounds(dist0, bwcheck)
			b1 = rd2d_kth_bounds(dist1, bwcheck)
			bwmin0 = b0[1]
			bwmax0 = b0[2]
			bwmin1 = b1[1]
			bwmax1 = b1[2]
			dn0 = min((max((dn0, bwmin0)), bwmax0))
			dn1 = min((max((dn1, bwmin1)), bwmax1))
		}
		if (covadj) {
			w0dn = rd2d_location_weights(centered, (dn0, dn0), kernel, kerneltype) :* side0
			w1dn = rd2d_location_weights(centered, (dn1, dn1), kernel, kerneltype) :* side1
		}
		else {
			w0dn = rd2d_location_weights(x0, (dn0, dn0), kernel, kerneltype)
			w1dn = rd2d_location_weights(x1, (dn1, dn1), kernel, kerneltype)
		}
		eNdn0 = sum(w0dn :> 0)
		eNdn1 = sum(w1dn :> 0)

		ybase = y
		fuzzybase = fuzzy
		if (covadj) {
			if (rows(fuzzy) == n & bwparam == "main") outcomes = y, fuzzy
			else outcomes = y
			design = rd2d_basis2(centered, p)
			w0full = w0dn
			w1full = w1dn
			gamma = rd2d_covariate_gamma(design, w0full, w1full, outcomes, covs)
			adjusted = rd2d_apply_covariates(outcomes, covs, gamma)
			ybase = adjusted[,1]
			if (rows(fuzzy) == n & bwparam == "main") fuzzybase = adjusted[,2]
		}

		ybw = ybase
		if (rows(fuzzy) == n & bwparam == "main") {
			outcomes = ybase, fuzzybase
			beta0 = rd2d_lm_multi_exact2(select(centered, side0), select(outcomes, side0),
				select(centered_dist, side0), dn0, p, kernel, kerneltype)
			beta1 = rd2d_lm_multi_exact2(select(centered, side1), select(outcomes, side1),
				select(centered_dist, side1), dn1, p, kernel, kerneltype)
			tauitt = (target * beta1[,1] - target * beta0[,1])[1,1]
			taufs = (target * beta1[,2] - target * beta0[,2])[1,1]
			if (taufs < . & abs(taufs) > sqrt(epsilon(1))) {
				graditt = 1 / taufs
				gradfs = -tauitt / taufs^2
				ybw = graditt :* ybase + gradfs :* fuzzybase
			}
		}
		ybw0 = select(ybw, side0)
		ybw1 = select(ybw, side1)

		vecq0 = rd2d_get_coeff_exact2(x0, dist0, target, p, dn0, kernel, kerneltype)
		vecq1 = rd2d_get_coeff_exact2(x1, dist1, target, p, dn1, kernel, kerneltype)
		thr0 = rd2d_quantile(dist0, .5)
		thr1 = rd2d_quantile(dist1, .5)
		if (covadj) {
			w0thr = rd2d_location_weights(centered, (thr0, thr0), kernel, kerneltype) :* side0
			w1thr = rd2d_location_weights(centered, (thr1, thr1), kernel, kerneltype) :* side1
		}
		else {
			w0thr = rd2d_location_weights(x0, (thr0, thr0), kernel, kerneltype)
			w1thr = rd2d_location_weights(x1, (thr1, thr1), kernel, kerneltype)
		}
		eNthr0 = sum(w0thr :> 0)
		eNthr1 = sum(w1thr :> 0)
		bn0 = thr0
		bn1 = thr1
		if (method == "dpi") {
			ybnv = ybw
			ybnb = ybw
			n_cov = 0
			if (covadj) {
				design = rd2d_basis2(centered, p + 1)
				w0full = w0dn
				w1full = w1dn
				gamma = rd2d_covariate_gamma(design, w0full, w1full, ybw, covs)
				ybnv = rd2d_apply_covariates(ybw, covs, gamma)
				n_cov = rd2d_covariate_rank(design, w0full, w1full, covs)
			}
			eN0 = eNdn0
			eN1 = eNdn1
			scales_v = rd2d_bw_vce_scales(fitmethod, vce, covadj, eN0, eN1, rd2d_basis_count2(p + 1), n_cov)
			n_cov = 0
			if (covadj) {
				design = rd2d_basis2(centered, p + 2)
				w0full = w0thr
				w1full = w1thr
				gamma = rd2d_covariate_gamma(design, w0full, w1full, ybw, covs)
				ybnb = rd2d_apply_covariates(ybw, covs, gamma)
				n_cov = rd2d_covariate_rank(design, w0full, w1full, covs)
			}
			eN0 = eNthr0
			eN1 = eNthr1
			scales_b = rd2d_bw_vce_scales(fitmethod, vce, covadj, eN0, eN1, rd2d_basis_count2(p + 2), n_cov)
			if (covadj) {
				ybnv0 = select(ybnv, side0)
				ybnv1 = select(ybnv, side1)
				ybnb0 = select(ybnb, side0)
				ybnb1 = select(ybnb, side1)
			}
			else {
				ybnv0 = ybw0
				ybnv1 = ybw1
				ybnb0 = ybw0
				ybnb1 = ybw1
			}
			bnconst0 = rd2d_bw_v2_exact2_adj(x0, ybw0, dist0,
				p + 1, vecq0, dn0, thr0, ., vce, kernel, kerneltype,
				cluster0,
				ybnv0, ybnb0, ybnb0,
				scales_v[1], scales_b[1])
			bnconst1 = rd2d_bw_v2_exact2_adj(x1, ybw1, dist1,
				p + 1, vecq1, dn1, thr1, ., vce, kernel, kerneltype,
				cluster1,
				ybnv1, ybnb1, ybnb1,
				scales_v[2], scales_b[2])
			bn0 = ((2 + 2*(p+1)) * bnconst0.V / (2 * (bnconst0.B^2 + scaleregul * bnconst0.Reg1)))^(1/(2*p + 6))
			bn1 = ((2 + 2*(p+1)) * bnconst1.V / (2 * (bnconst1.B^2 + scaleregul * bnconst1.Reg1)))^(1/(2*p + 6))
			if (bwcheck > 0) {
				bn0 = min((max((bn0, bwmin0)), bwmax0))
				bn1 = min((max((bn1, bwmin1)), bwmax1))
			}
		}

		yhnv = ybw
		yhnb = ybw
		yhnt = ybw
		n_cov = 0
		if (covadj) {
			design = rd2d_basis2(centered, p)
			w0full = w0dn
			w1full = w1dn
			gamma = rd2d_covariate_gamma(design, w0full, w1full, ybw, covs)
			yhnv = rd2d_apply_covariates(ybw, covs, gamma)
			n_cov = rd2d_covariate_rank(design, w0full, w1full, covs)
		}
		eN0 = eNdn0
		eN1 = eNdn1
		scales_v = rd2d_bw_vce_scales(fitmethod, vce, covadj, eN0, eN1, rd2d_basis_count2(p), n_cov)
		if (covadj) {
			w0bn = rd2d_location_weights(centered, (bn0, bn0), kernel, kerneltype) :* side0
			w1bn = rd2d_location_weights(centered, (bn1, bn1), kernel, kerneltype) :* side1
		}
		else {
			w0bn = rd2d_location_weights(x0, (bn0, bn0), kernel, kerneltype)
			w1bn = rd2d_location_weights(x1, (bn1, bn1), kernel, kerneltype)
		}
		eNbn0 = sum(w0bn :> 0)
		eNbn1 = sum(w1bn :> 0)
		n_cov = 0
		if (covadj) {
			design = rd2d_basis2(centered, p + 1)
			w0full = w0bn
			w1full = w1bn
			gamma = rd2d_covariate_gamma(design, w0full, w1full, ybw, covs)
			yhnb = rd2d_apply_covariates(ybw, covs, gamma)
			n_cov = rd2d_covariate_rank(design, w0full, w1full, covs)
		}
		eN0 = eNbn0
		eN1 = eNbn1
		scales_b = rd2d_bw_vce_scales(fitmethod, vce, covadj, eN0, eN1, rd2d_basis_count2(p + 1), n_cov)
		if (covadj) {
			design = rd2d_basis2(centered, p + 2)
			w0full = w0thr
			w1full = w1thr
			gamma = rd2d_covariate_gamma(design, w0full, w1full, ybw, covs)
			yhnt = rd2d_apply_covariates(ybw, covs, gamma)
		}
		if (covadj) {
			yhnv0 = select(yhnv, side0)
			yhnv1 = select(yhnv, side1)
			yhnb0 = select(yhnb, side0)
			yhnb1 = select(yhnb, side1)
			yhnt0 = select(yhnt, side0)
			yhnt1 = select(yhnt, side1)
		}
		else {
			yhnv0 = ybw0
			yhnv1 = ybw1
			yhnb0 = ybw0
			yhnb1 = ybw1
			yhnt0 = ybw0
			yhnt1 = ybw1
		}
		hnconst0 = rd2d_bw_v2_exact2_adj(x0, ybw0, dist0,
			p, target, dn0, bn0, thr0, vce, kernel, kerneltype,
			cluster0,
			yhnv0, yhnb0, yhnt0,
			scales_v[1], scales_b[1])
		hnconst1 = rd2d_bw_v2_exact2_adj(x1, ybw1, dist1,
			p, target, dn1, bn1, thr1, vce, kernel, kerneltype,
			cluster1,
			yhnv1, yhnb1, yhnt1,
			scales_v[2], scales_b[2])

		if (rd2d_bwbase(bwselect) == "mserd" | rd2d_bwbase(bwselect) == "imserd") {
			hn = ((2 + 2*derivsum) * (hnconst0.V + hnconst1.V) /
				(derivdenom * ((hnconst0.B + scalebiascrct*hnconst0.Reg2 -
				hnconst1.B - scalebiascrct*hnconst1.Reg2)^2 +
				scaleregul*hnconst0.Reg1 + scaleregul*hnconst1.Reg1)))^(1/(2*p + 4))
			if (bwcheck > 0) hn = min((max((hn, bwmin0, bwmin1)), max((bwmax0, bwmax1))))
			hn0 = hn
			hn1 = hn
		}
		else {
			hn0 = ((2 + 2*derivsum) * hnconst0.V /
				(derivdenom * ((hnconst0.B + scalebiascrct*hnconst0.Reg2)^2 +
				scaleregul*hnconst0.Reg1)))^(1/(2*p + 4))
			hn1 = ((2 + 2*derivsum) * hnconst1.V /
				(derivdenom * ((hnconst1.B + scalebiascrct*hnconst1.Reg2)^2 +
				scaleregul*hnconst1.Reg1)))^(1/(2*p + 4))
			if (bwcheck > 0) {
				hn0 = min((max((hn0, bwmin0)), bwmax0))
				hn1 = min((max((hn1, bwmin1)), bwmax1))
			}
		}

		results[j,.] = (bwork[j,1], bwork[j,2], hn0, hn0, hn1, hn1,
			eNdn0, eNdn1,
			hnconst0.B, hnconst1.B, hnconst0.V, hnconst1.V,
			hnconst0.Reg2, hnconst1.Reg2, hnconst0.Reg1, hnconst1.Reg1)
		bw[j,.] = (b[j,1], b[j,2], hn0*sd[1], hn0*sd[2],
			hn1*sd[1], hn1*sd[2], results[j,7], results[j,8])
	}

	if (rd2d_bwbase(bwselect) == "imserd") {
		VV = mean(results[,11]) + mean(results[,12])
		BB = mean((results[,9] :+ scalebiascrct:*results[,13] :-
			results[,10] :- scalebiascrct:*results[,14]):^2 :+
			scaleregul:*results[,15] :+ scaleregul:*results[,16])
		hn = ((2 + 2*derivsum) * VV / (derivdenom * BB))^(1/(2*p + 4))
		for (i=1; i<=neval; i++) bw[i,3..6] = (hn*sd[1], hn*sd[2], hn*sd[1], hn*sd[2])
	}
	else if (rd2d_bwbase(bwselect) == "imsetwo") {
		VV = mean(results[,11])
		BB0 = mean((results[,9] :+ scalebiascrct:*results[,13]):^2 :+
			scaleregul:*results[,15])
		hn0 = ((2 + 2*derivsum) * VV / (derivdenom * BB0))^(1/(2*p + 4))
		VV = mean(results[,12])
		BB1 = mean((results[,10] :+ scalebiascrct:*results[,14]):^2 :+
			scaleregul:*results[,16])
		hn1 = ((2 + 2*derivsum) * VV / (derivdenom * BB1))^(1/(2*p + 4))
		for (i=1; i<=neval; i++) bw[i,3..6] = (hn0*sd[1], hn0*sd[2], hn1*sd[1], hn1*sd[2])
	}

	if (rd2d_is_cer(bwselect)) {
		if (rd2d_is_common(bwselect)) {
			factor0 = rd2d_cer_factor(M, p)
			bw[,3..6] = bw[,3..6] :* factor0
		}
		else {
			factor0 = rd2d_cer_factor(M0, p)
			factor1 = rd2d_cer_factor(M1, p)
			bw[,3..4] = bw[,3..4] :* factor0
			bw[,5..6] = bw[,5..6] :* factor1
		}
	}
	return(bw)
}

struct rd2d_orderfit scalar rd2d_fit_location_order(real matrix x, real colvector d,
	real matrix b, real matrix outcomes, real matrix H, real scalar p,
	real rowvector deriv, real matrix tangvec, string scalar kernel,
	string scalar kerneltype, string scalar vce, real colvector cluster,
	string scalar fitmethod, real matrix covs)
{
	struct rd2d_orderfit scalar out
	struct rd2d_lfit scalar fit0, fit1
	real matrix centered, centered0, centered1, design, design0, design1
	real matrix outcomes_fit, outcomes0, outcomes1, gamma
	real colvector w0, w1, active, clactive, cl0, cl1, cluster_groups
	real colvector side0, side1, idx0, idx1
	real rowvector target
	real scalar j, n, neval, nout, ninf, covadj, n_cov, clustered, sidefit, e0, e1, kpoly
	real scalar scale0, scale1
	string scalar vce_local

	n = rows(x)
	neval = rows(b)
	nout = cols(outcomes)
	covadj = (cols(covs) > 0)
	clustered = (rows(cluster) == n)
	sidefit = (!covadj & !clustered & neval <= 6)
	cluster_groups = clustered ? uniqrows(sort(cluster,1)) : J(0,1,.)
	ninf = (clustered ? rows(cluster_groups) : n)
	fitmethod = strlower(strtrim(fitmethod))
	vce = strlower(strtrim(vce))
	out.mu0 = J(neval, nout, .)
	out.mu1 = J(neval, nout, .)
	out.se0 = J(neval, nout, .)
	out.se1 = J(neval, nout, .)
	out.N0 = J(neval, 1, 0)
	out.N1 = J(neval, 1, 0)
	out.infl0_1 = J(neval, ninf, 0)
	out.infl1_1 = J(neval, ninf, 0)
	out.infl0_2 = J(neval, ninf, 0)
	out.infl1_2 = J(neval, ninf, 0)
	if (sidefit) {
		side0 = d :== 0
		side1 = d :!= 0
		idx0 = selectindex(side0)
		idx1 = selectindex(side1)
		outcomes0 = outcomes[idx0,.]
		outcomes1 = outcomes[idx1,.]
	}

	for (j=1; j<=neval; j++) {
		target = rd2d_target2(p, deriv, tangvec, j)
		n_cov = 0
		if (sidefit) {
			centered0 = x[idx0,.] :- (J(rows(idx0),1,1)*b[j,.])
			centered1 = x[idx1,.] :- (J(rows(idx1),1,1)*b[j,.])
			design0 = rd2d_basis2(centered0, p)
			design1 = rd2d_basis2(centered1, p)
			w0 = rd2d_location_weights(centered0, H[j,(1,2)], kernel, kerneltype)
			w1 = rd2d_location_weights(centered1, H[j,(3,4)], kernel, kerneltype)
			e0 = sum(w0 :> 0)
			e1 = sum(w1 :> 0)
			kpoly = cols(design0)
			vce_local = ((fitmethod == "joint" | covadj) & vce == "hc1") ? "hc0" : vce
			scale0 = .
			scale1 = .
			if ((fitmethod == "joint" | covadj) & (vce == "hc1" | clustered)) {
				if (fitmethod == "joint") {
					scale0 = rd2d_joint_vce_scale(vce, e0 + e1, 2*kpoly + n_cov, J(0,1,.), clustered)
					scale1 = scale0
				}
				else {
					scale0 = rd2d_joint_vce_scale(vce, e0, kpoly + n_cov, J(0,1,.), clustered)
					scale1 = rd2d_joint_vce_scale(vce, e1, kpoly + n_cov, J(0,1,.), clustered)
				}
			}
			fit0 = rd2d_local_fit_scaled(design0, w0, outcomes0, target, vce_local, J(0,1,.), scale0, cluster_groups)
			fit1 = rd2d_local_fit_scaled(design1, w1, outcomes1, target, vce_local, J(0,1,.), scale1, cluster_groups)
		}
		else {
			centered = x :- (J(n,1,1)*b[j,.])
			design = rd2d_basis2(centered, p)
			w0 = rd2d_location_weights(centered, H[j,(1,2)], kernel, kerneltype) :* (d:==0)
			w1 = rd2d_location_weights(centered, H[j,(3,4)], kernel, kerneltype) :* (d:!=0)
			outcomes_fit = outcomes
			if (covadj) {
				gamma = rd2d_covariate_gamma(design, w0, w1, outcomes, covs)
				n_cov = rd2d_covariate_rank(design, w0, w1, covs)
				outcomes_fit = rd2d_apply_covariates(outcomes, covs, gamma)
			}
			e0 = sum(w0 :> 0)
			e1 = sum(w1 :> 0)
			kpoly = cols(design)
			vce_local = ((fitmethod == "joint" | covadj) & vce == "hc1") ? "hc0" : vce
			scale0 = .
			scale1 = .
			if ((fitmethod == "joint" | covadj) & (vce == "hc1" | clustered)) {
				if (fitmethod == "joint") {
					active = (w0 :> 0) :| (w1 :> 0)
					clactive = clustered ? select(cluster, active) : J(0,1,.)
					scale0 = rd2d_joint_vce_scale(vce, e0 + e1, 2*kpoly + n_cov, clactive, clustered)
					scale1 = scale0
				}
				else {
					cl0 = clustered ? select(cluster, w0 :> 0) : J(0,1,.)
					cl1 = clustered ? select(cluster, w1 :> 0) : J(0,1,.)
					scale0 = rd2d_joint_vce_scale(vce, e0, kpoly + n_cov, cl0, clustered)
					scale1 = rd2d_joint_vce_scale(vce, e1, kpoly + n_cov, cl1, clustered)
				}
			}
			fit0 = rd2d_local_fit_scaled(design, w0, outcomes_fit, target, vce_local, cluster, scale0, cluster_groups)
			fit1 = rd2d_local_fit_scaled(design, w1, outcomes_fit, target, vce_local, cluster, scale1, cluster_groups)
		}
		out.mu0[j,.] = fit0.estimate
		out.mu1[j,.] = fit1.estimate
		out.se0[j,.] = fit0.se
		out.se1[j,.] = fit1.se
		out.N0[j] = fit0.n_eff
		out.N1[j] = fit1.n_eff
		if (sidefit) {
			out.infl0_1[j,idx0] = fit0.influence[,1]'
			out.infl1_1[j,idx1] = fit1.influence[,1]'
			if (nout >= 2) {
				out.infl0_2[j,idx0] = fit0.influence[,2]'
				out.infl1_2[j,idx1] = fit1.influence[,2]'
			}
		}
		else {
			out.infl0_1[j,.] = fit0.influence[,1]'
			out.infl1_1[j,.] = fit1.influence[,1]'
			if (nout >= 2) {
				out.infl0_2[j,.] = fit0.influence[,2]'
				out.infl1_2[j,.] = fit1.influence[,2]'
			}
		}
	}
	return(out)
}

real matrix rd2d_ci_table(real matrix b, real colvector estp, real colvector sep,
	real colvector estq, real colvector seq, real matrix H, real colvector N0,
	real colvector N1, real scalar level, string scalar side, real scalar distance)
{
	real matrix out
	real colvector tval, pval, cil, ciu
	real scalar cval

	tval = estq :/ seq
	pval = 2 :* (1 :- normal(abs(tval)))
	side = strlower(strtrim(side))
	if (side == "left" | side == "right") cval = invnormal(level/100)
	else cval = invnormal((level+100)/200)
	if (side == "left") {
		cil = J(rows(estq),1,.)
		ciu = estq :+ cval:*seq
	}
	else if (side == "right") {
		cil = estq :- cval:*seq
		ciu = J(rows(estq),1,.)
	}
	else {
		cil = estq :- cval:*seq
		ciu = estq :+ cval:*seq
	}
	out = b, estp, sep, estq, seq, tval, pval, cil, ciu, H, N0, N1
	return(out)
}

void rd2d_location_stata(string scalar main_name, string scalar bw_name,
	string scalar main0_name, string scalar main1_name, string scalar itt_name,
	string scalar fs_name, string scalar itt0_name, string scalar itt1_name,
	string scalar fs0_name, string scalar fs1_name, string scalar vmain_name,
	string scalar vitt_name, string scalar vfs_name, string scalar vitt0_name,
	string scalar vitt1_name, string scalar vfs0_name, string scalar vfs1_name,
	string scalar yvar,
	string scalar x1var, string scalar x2var, string scalar dvar,
	string scalar fuzzyvar, string scalar clustvar, string scalar touse,
	string scalar bstr, string scalar hstr, string scalar derivstr,
	string scalar tangstr, real scalar p, real scalar q, string scalar kernel,
	string scalar kerneltype, string scalar vce, real scalar level,
	string scalar side, string scalar bwselect, string scalar bwparam,
	string scalar method, real scalar bwcheck, string scalar masspoints,
	real scalar scaleregul, real scalar scalebiascrct, real scalar stdvars,
	string scalar fitmethod, string scalar covsvars)
{
	real matrix data, x, b, bw, H, outcomes, covs
	real matrix main, main0, main1, itt, fs, itt0, itt1, fs0, fs1
	real matrix Vmain, Vitt, Vfs, Vitt0, Vitt1, Vfs0, Vfs1
	real colvector y, d, fuzzy, cluster, estp, estq, sep, seq, taup, tauq
	real colvector itt_p, itt_q, fs_p, fs_q, se_itt_p, se_itt_q, se_fs_p, se_fs_q
	real matrix inflp, inflq, infl_itt_p, infl_itt_q, infl_fs_p, infl_fs_q
	real matrix deriv_raw
	real rowvector deriv
	real scalar covstart, covend
	real matrix tangvec
	struct rd2d_orderfit scalar fitp, fitq

	if (q < 0) q = p + 1
	b = rd2d_parse_numlist(bstr, 2)
	deriv_raw = rd2d_parse_numlist(derivstr, 1)
	deriv = rows(deriv_raw) == 0 ? (0,0) : vec(deriv_raw)'
	if (cols(deriv) == 1) deriv = (deriv[1], 0)
	if (cols(deriv) != 2) {
		errprintf("deriv() must contain one or two nonnegative integers\n")
		_error(198)
	}
	tangvec = rd2d_parse_numlist(tangstr, 2)

	data = st_data(., yvar+" "+x1var+" "+x2var+" "+dvar+
		(fuzzyvar=="" ? "" : " "+fuzzyvar)+
		(covsvars=="" ? "" : " "+covsvars)+
		(clustvar=="" ? "" : " "+clustvar), touse)
	y = data[,1]
	x = data[,(2,3)]
	d = data[,4]
	if (fuzzyvar != "") fuzzy = data[,5]
	else fuzzy = J(0,1,.)
	if (covsvars != "") {
		covstart = 5 + (fuzzyvar != "")
		covend = covstart + cols(tokens(covsvars)) - 1
		covs = data[,covstart..covend]
	}
	else covs = J(rows(data), 0, .)
	if (clustvar != "") cluster = data[,cols(data)]
	else cluster = J(0,1,.)

	fitmethod = (strlower(strtrim(fitmethod)) == "separate" ? "separate" : "joint")
	bw = rd2d_location_bw_full(y, x, d, b, hstr, fuzzy, cluster, p,
		deriv, tangvec, kernel, kerneltype, bwselect, bwparam, method,
		vce, bwcheck, masspoints, scaleregul, scalebiascrct, stdvars,
		fitmethod, covs)
	H = bw[,(3,4,5,6)]
	H = rd2d_apply_location_bwcheck(H, x, d, b, bwcheck, kerneltype)
	if (fuzzyvar == "") outcomes = y
	else outcomes = y, fuzzy
	fitp = rd2d_fit_location_order(x, d, b, outcomes, H, p, deriv, tangvec, kernel, kerneltype, vce, cluster, fitmethod, covs)
	fitq = (q == p ? fitp : rd2d_fit_location_order(x, d, b, outcomes, H, q, deriv, tangvec, kernel, kerneltype, vce, cluster, fitmethod, covs))

	if (fuzzyvar == "") {
		taup = fitp.mu1[,1] - fitp.mu0[,1]
		tauq = fitq.mu1[,1] - fitq.mu0[,1]
		if (rows(cluster) == rows(y) & fitmethod == "separate") {
			inflp = fitp.infl1_1, -fitp.infl0_1
			inflq = fitq.infl1_1, -fitq.infl0_1
		}
		else {
			inflp = fitp.infl1_1 - fitp.infl0_1
			inflq = fitq.infl1_1 - fitq.infl0_1
		}
		sep = sqrt(diagonal(inflp*inflp'))
		seq = sqrt(diagonal(inflq*inflq'))
		main = rd2d_ci_table(b, taup, sep, tauq, seq, H, fitp.N0, fitp.N1, level, side, 0)
		main0 = rd2d_ci_table(b, fitp.mu0[,1], fitp.se0[,1], fitq.mu0[,1], fitq.se0[,1],
			(H[,(1,2)], J(rows(H),2,.)), fitp.N0, J(rows(H),1,.), level, side, 0)
		main1 = rd2d_ci_table(b, fitp.mu1[,1], fitp.se1[,1], fitq.mu1[,1], fitq.se1[,1],
			(J(rows(H),2,.), H[,(3,4)]), J(rows(H),1,.), fitp.N1, level, side, 0)
		itt = J(0,0,.)
		fs = J(0,0,.)
		itt0 = J(0,0,.)
		itt1 = J(0,0,.)
		fs0 = J(0,0,.)
		fs1 = J(0,0,.)
		Vmain = inflq * inflq'
		Vitt = J(0,0,.)
		Vfs = J(0,0,.)
		Vitt0 = J(0,0,.)
		Vitt1 = J(0,0,.)
		Vfs0 = J(0,0,.)
		Vfs1 = J(0,0,.)
	}
	else {
		itt_p = fitp.mu1[,1] - fitp.mu0[,1]
		itt_q = fitq.mu1[,1] - fitq.mu0[,1]
		fs_p = fitp.mu1[,2] - fitp.mu0[,2]
		fs_q = fitq.mu1[,2] - fitq.mu0[,2]
		taup = itt_p :/ fs_p
		tauq = itt_q :/ fs_q
		if (rows(cluster) == rows(y) & fitmethod == "separate") {
			infl_itt_p = fitp.infl1_1, -fitp.infl0_1
			infl_itt_q = fitq.infl1_1, -fitq.infl0_1
			infl_fs_p = fitp.infl1_2, -fitp.infl0_2
			infl_fs_q = fitq.infl1_2, -fitq.infl0_2
		}
		else {
			infl_itt_p = fitp.infl1_1 - fitp.infl0_1
			infl_itt_q = fitq.infl1_1 - fitq.infl0_1
			infl_fs_p = fitp.infl1_2 - fitp.infl0_2
			infl_fs_q = fitq.infl1_2 - fitq.infl0_2
		}
		inflp = (infl_itt_p :/ (fs_p*J(1,cols(infl_itt_p),1))) :-
			((itt_p:/fs_p:^2)*J(1,cols(infl_fs_p),1)) :* infl_fs_p
		inflq = (infl_itt_q :/ (fs_q*J(1,cols(infl_itt_q),1))) :-
			((itt_q:/fs_q:^2)*J(1,cols(infl_fs_q),1)) :* infl_fs_q
		sep = sqrt(diagonal(inflp*inflp'))
		seq = sqrt(diagonal(inflq*inflq'))
		se_itt_p = sqrt(diagonal(infl_itt_p*infl_itt_p'))
		se_itt_q = sqrt(diagonal(infl_itt_q*infl_itt_q'))
		se_fs_p = sqrt(diagonal(infl_fs_p*infl_fs_p'))
		se_fs_q = sqrt(diagonal(infl_fs_q*infl_fs_q'))
		main = rd2d_ci_table(b, taup, sep, tauq, seq, H, fitp.N0, fitp.N1, level, side, 0)
		itt = rd2d_ci_table(b, itt_p, se_itt_p, itt_q, se_itt_q, H, fitp.N0, fitp.N1, level, side, 0)
		fs = rd2d_ci_table(b, fs_p, se_fs_p, fs_q, se_fs_q, H, fitp.N0, fitp.N1, level, side, 0)
		itt0 = rd2d_ci_table(b, fitp.mu0[,1], fitp.se0[,1], fitq.mu0[,1], fitq.se0[,1],
			(H[,(1,2)], J(rows(H),2,.)), fitp.N0, J(rows(H),1,.), level, side, 0)
		itt1 = rd2d_ci_table(b, fitp.mu1[,1], fitp.se1[,1], fitq.mu1[,1], fitq.se1[,1],
			(J(rows(H),2,.), H[,(3,4)]), J(rows(H),1,.), fitp.N1, level, side, 0)
		fs0 = rd2d_ci_table(b, fitp.mu0[,2], fitp.se0[,2], fitq.mu0[,2], fitq.se0[,2],
			(H[,(1,2)], J(rows(H),2,.)), fitp.N0, J(rows(H),1,.), level, side, 0)
		fs1 = rd2d_ci_table(b, fitp.mu1[,2], fitp.se1[,2], fitq.mu1[,2], fitq.se1[,2],
			(J(rows(H),2,.), H[,(3,4)]), J(rows(H),1,.), fitp.N1, level, side, 0)
		main0 = J(0,0,.)
		main1 = J(0,0,.)
		Vmain = inflq * inflq'
		Vitt = infl_itt_q * infl_itt_q'
		Vfs = infl_fs_q * infl_fs_q'
		Vitt0 = fitq.infl0_1 * fitq.infl0_1'
		Vitt1 = fitq.infl1_1 * fitq.infl1_1'
		Vfs0 = fitq.infl0_2 * fitq.infl0_2'
		Vfs1 = fitq.infl1_2 * fitq.infl1_2'
	}
	bw = b, H, fitp.N0, fitp.N1
	st_matrix(main_name, main)
	st_matrix(bw_name, bw)
	st_matrix(main0_name, main0)
	st_matrix(main1_name, main1)
	st_matrix(itt_name, itt)
	st_matrix(fs_name, fs)
	st_matrix(itt0_name, itt0)
	st_matrix(itt1_name, itt1)
	st_matrix(fs0_name, fs0)
	st_matrix(fs1_name, fs1)
	st_matrix(vmain_name, Vmain)
	st_matrix(vitt_name, Vitt)
	st_matrix(vfs_name, Vfs)
	st_matrix(vitt0_name, Vitt0)
	st_matrix(vitt1_name, Vitt1)
	st_matrix(vfs0_name, Vfs0)
	st_matrix(vfs1_name, Vfs1)
}

void rdbw2d_location_stata(string scalar bw_name, string scalar yvar,
	string scalar x1var, string scalar x2var, string scalar dvar,
	string scalar fuzzyvar, string scalar clustvar, string scalar touse,
	string scalar bstr, string scalar derivstr, string scalar tangstr,
	real scalar p, string scalar kernel,
	string scalar kerneltype, string scalar bwselect, string scalar bwparam,
	string scalar method, string scalar vce, real scalar bwcheck,
	string scalar masspoints, real scalar scaleregul, real scalar scalebiascrct,
	real scalar stdvars, string scalar fitmethod, string scalar covsvars)
{
	real matrix data, x, b, bw, tangvec, deriv_raw, covs
	real colvector y, d, fuzzy, cluster
	real rowvector deriv
	real scalar covstart, covend

	b = rd2d_parse_numlist(bstr, 2)
	deriv_raw = rd2d_parse_numlist(derivstr, 1)
	deriv = rows(deriv_raw) == 0 ? (0,0) : vec(deriv_raw)'
	if (cols(deriv) == 1) deriv = (deriv[1], 0)
	tangvec = rd2d_parse_numlist(tangstr, 2)
	data = st_data(., yvar+" "+x1var+" "+x2var+" "+dvar+
		(fuzzyvar=="" ? "" : " "+fuzzyvar)+
		(covsvars=="" ? "" : " "+covsvars)+
		(clustvar=="" ? "" : " "+clustvar), touse)
	y = data[,1]
	x = data[,(2,3)]
	d = data[,4]
	if (fuzzyvar != "") fuzzy = data[,5]
	else fuzzy = J(0,1,.)
	if (covsvars != "") {
		covstart = 5 + (fuzzyvar != "")
		covend = covstart + cols(tokens(covsvars)) - 1
		covs = data[,covstart..covend]
	}
	else covs = J(rows(data), 0, .)
	if (clustvar != "") cluster = data[,cols(data)]
	else cluster = J(0,1,.)
	fitmethod = (strlower(strtrim(fitmethod)) == "separate" ? "separate" : "joint")
	bw = rd2d_location_bw_full(y, x, d, b, "", fuzzy, cluster, p,
		deriv, tangvec, kernel, kerneltype, bwselect, bwparam, method,
		vce, bwcheck, masspoints, scaleregul, scalebiascrct, stdvars,
		fitmethod, covs)
	st_matrix(bw_name, bw)
}

real matrix rd2d_distance_bw_simple(real matrix dist, real matrix b, string scalar hstr,
	real scalar p, string scalar kernel, string scalar bwselect, real scalar bwcheck,
	real scalar kink1)
{
	real matrix H, hraw, bw
	real colvector sdist, u, u0, u1
	real rowvector b0, b1
	real scalar j, neval, hrot, factor, M0, M1, M, user_h

	neval = cols(dist)
	hraw = rd2d_parse_numlist(hstr, 1)
	user_h = rows(hraw) > 0
	if (rows(hraw) == 1 & cols(hraw) == 1) H = J(neval, 2, hraw[1,1])
	else if (rows(hraw) > 0) {
		if (rows(hraw)*cols(hraw) != neval*2) {
			errprintf("h() must be a scalar or contain 2 values per evaluation point\n")
			_error(198)
		}
		H = rd2d_parse_numlist(hstr, 2)
	}
	else {
		H = J(neval, 2, .)
		for (j=1; j<=neval; j++) {
			sdist = dist[,j]
			u = abs(sdist)
			hrot = rd2d_rot_distance(u, kernel)
			H[j,.] = (hrot, hrot)
		}
		if (rd2d_is_cer(bwselect)) {
			for (j=1; j<=neval; j++) {
				sdist = dist[,j]
				M0 = rows(uniqrows(select(sdist, sdist:<0)))
				M1 = rows(uniqrows(select(sdist, sdist:>=0)))
				M = M0 + M1
				if (rd2d_is_common(bwselect)) H[j,.] = H[j,.] :* rd2d_cer_factor(M, p)
				else {
					H[j,1] = H[j,1] * rd2d_cer_factor(M0, p)
					H[j,2] = H[j,2] * rd2d_cer_factor(M1, p)
				}
			}
		}
	}
	for (j=1; j<=neval; j++) {
		sdist = dist[,j]
		u0 = sort(abs(select(sdist, sdist:<0)), 1)
		u1 = sort(abs(select(sdist, sdist:>=0)), 1)
		if (bwcheck > 0) {
			b0 = rd2d_kth_bounds(u0, bwcheck)
			b1 = rd2d_kth_bounds(u1, bwcheck)
			H[j,1] = min((max((H[j,1], b0[1])), b0[2]))
			H[j,2] = min((max((H[j,2], b1[1])), b1[2]))
		}
		if (kink1 & !user_h) {
			M0 = max((rows(u0),1))
			M1 = max((rows(u1),1))
			if (rd2d_is_common(bwselect)) {
				factor = rows(dist)^(1/(2*p+4)-1/4)
				H[j,.] = H[j,.] :* factor
			}
			else {
				H[j,1] = H[j,1] * M0^(1/(2*p+4)-1/4)
				H[j,2] = H[j,2] * M1^(1/(2*p+4)-1/4)
			}
		}
	}
	bw = b, H, J(neval,1,0), J(neval,1,0)
	for (j=1; j<=neval; j++) {
		sdist = dist[,j]
		bw[j,5] = sum((sdist:<0) :& (rd2d_distance_weights(sdist, H[j,1], 0, kernel):>0))
		bw[j,6] = sum((sdist:>=0) :& (rd2d_distance_weights(sdist, H[j,2], 1, kernel):>0))
	}
	return(bw)
}

struct rd2d_polyfit scalar rd2d_poly_fit1(real colvector y, real colvector u,
	real scalar degree)
{
	struct rd2d_polyfit scalar out
	real matrix X, invXX
	real colvector res
	real scalar sigma2

	X = rd2d_basis1(u, degree)
	invXX = pinv(X' * X)
	out.coef = invXX * (X' * y)
	res = y - X * out.coef
	sigma2 = quadcross(res, res) / max((rows(X) - cols(X), 1))
	out.se = sqrt(diagonal(invXX) :* sigma2)
	return(out)
}

real matrix rd2d_vce_const1(real matrix wR, real colvector resd, real scalar h)
{
	real matrix scores
	scores = wR :* (resd * J(1, cols(wR), 1))
	return(scores' * scores * h^2)
}

real rowvector rd2d_distance_intercepts_multi(real colvector y, real colvector fuzzy,
	real colvector u, real scalar h, real scalar p, string scalar kernel)
{
	real colvector w, idx, ew, sqrtw
	real matrix R, sqrtw_R, invG, outcomes, beta

	w = rd2d_kernel(u:/h, kernel) :/ max((h^2, epsilon(1)))
	idx = selectindex(w :> 0)
	if (rows(idx) <= p + 1) return((., .))
	ew = w[idx]
	R = rd2d_basis1(u[idx]:/h, p)
	sqrtw = sqrt(ew)
	sqrtw_R = R :* (sqrtw * J(1, cols(R), 1))
	invG = pinv(sqrtw_R' * sqrtw_R)
	outcomes = y[idx], fuzzy[idx]
	beta = invG * (R' * (outcomes :* (ew * J(1, cols(outcomes), 1))))
	return(beta[1,1], beta[1,2])
}

real rowvector rd2d_distance_side_constants_adj(real colvector u, real colvector y,
	real scalar h, real scalar p, struct rd2d_polyfit scalar model,
	string scalar vce, string scalar kernel, real scalar scale_override)
{
	real colvector w, idx, ew, eX, eY, resd, hii, pow
	real matrix R, sqrtw_R, invG, sigmahalf, sigma, pmatrix
	real rowvector vec, coeff
	real scalar eN, k, B, V, Reg

	w = rd2d_kernel(u:/h, kernel) :/ max((h^2, epsilon(1)))
	idx = selectindex(w :> 0)
	if (rows(idx) == 0) return((., ., ., 0))
	ew = w[idx]
	eX = u[idx]
	eY = y[idx]
	R = rd2d_basis1(eX:/h, p)
	sqrtw_R = R :* (sqrt(ew) * J(1, cols(R), 1))
	invG = pinv(sqrtw_R' * sqrtw_R)
	sigmahalf = R :* (ew * J(1, cols(R), 1))
	resd = eY - rd2d_basis1(eX, rows(model.coef)-1) * model.coef
	eN = rows(idx)
	k = p + 1
	vce = strlower(strtrim(vce))
	if (scale_override >= . & vce == "hc1" & eN > k) {
		resd = resd :* sqrt(eN / (eN - k))
	}
	else if (vce == "hc2" | vce == "hc3") {
		hii = rowsum((sqrtw_R * invG) :* sqrtw_R)
		if (vce == "hc2") resd = resd :* sqrt(1 :/ rowmax((1 :- hii, J(eN,1,1e-12))))
		else              resd = resd :* (1 :/ rowmax((1 :- hii, J(eN,1,1e-12))))
	}
	sigma = rd2d_vce_const1(sigmahalf, resd, h)
	vec = J(1, p + 1, 0)
	vec[1] = 1
	pow = (eX:/h):^(p + 1) :* ew
	pmatrix = R' * pow
	coeff = J(1, p + 2, 0)
	coeff[p + 2] = (vec * invG * pmatrix)[1,1]
	B = (coeff * model.coef)[1,1]
	Reg = coeff[p + 2]^2 * model.se[p + 2]^2
	V = (vec * invG * sigma * invG * vec')[1,1]
	if (scale_override < .) V = V * scale_override
	return((B, V, Reg, eN))
}

real rowvector rd2d_distance_side_constants(real colvector u, real colvector y,
	real scalar h, real scalar p, struct rd2d_polyfit scalar model,
	string scalar vce, string scalar kernel)
{
	return(rd2d_distance_side_constants_adj(u, y, h, p, model, vce, kernel, .))
}

real matrix rd2d_distance_bw_full(real matrix dist, real matrix b, string scalar hstr,
	real colvector y, real colvector fuzzy, real scalar p, string scalar kernel,
	string scalar bwselect, string scalar bwparam, string scalar vce,
	real scalar bwcheck, real scalar kink1, string scalar kinkposstr,
	real scalar scaleregul, real scalar cqt, string scalar fitmethod,
	real matrix covs)
{
	real matrix hraw, H, bw, results
	real colvector sdist, u, u0, u1, y0, y1, f0, f1, idx0, idx1
	real colvector ybase, fuzzybase
	real colvector kinkpos, distkink
	real matrix design, gamma, adjusted, outcomes
	real colvector w0full, w1full
	real rowvector b0, b1, c0, c1, mu0, mu1
	real rowvector scales
	struct rd2d_polyfit scalar fit0, fit1
	real scalar neval, j, dn, dn0, dn1, bwmin0, bwmin1, bwmax0, bwmax1
	real scalar thr0, thr1, tauitt, taufs, graditt, gradfs
	real scalar VV, BB, BB0, BB1, hn, h0, h1, smooth, factor, M0, M1, M
	real scalar covadj, n_cov, eN0, eN1

	hraw = rd2d_parse_numlist(hstr, 1)
	if (rows(hraw) > 0) {
		if (strtrim(kinkposstr) != "") {
			errprintf("kinkposition() applies only to automatic bandwidth selection; omit h() to use it\n")
			_error(198)
		}
		return(rd2d_distance_bw_simple(dist, b, hstr, p, kernel, bwselect,
			bwcheck, kink1))
	}

	neval = cols(dist)
	kinkpos = rd2d_parse_kink_position(kinkposstr, neval, b)
	if (sum(kinkpos) > 0 & kink1) {
		errprintf("use either kinkposition() or kinkunknown(), not both\n")
		_error(198)
	}
	results = J(neval, 14, .)
	bwselect = strlower(strtrim(bwselect))
	bwparam = strlower(strtrim(bwparam))
	fitmethod = (strlower(strtrim(fitmethod)) == "separate" ? "separate" : "joint")
	covadj = (cols(covs) > 0)

	for (j=1; j<=neval; j++) {
		sdist = dist[,j]
		u = abs(sdist)
		idx0 = sdist :< 0
		idx1 = sdist :>= 0
		u0 = select(u, idx0)
		u1 = select(u, idx1)
		ybase = y
		fuzzybase = fuzzy
		dn = rd2d_rot_distance(u, kernel)
		dn0 = dn
		dn1 = dn
		bwmin0 = 0
		bwmin1 = 0
		bwmax0 = max(u0)
		bwmax1 = max(u1)
		if (bwcheck > 0) {
			b0 = rd2d_kth_bounds(u0, bwcheck)
			b1 = rd2d_kth_bounds(u1, bwcheck)
			bwmin0 = b0[1]
			bwmax0 = b0[2]
			bwmin1 = b1[1]
			bwmax1 = b1[2]
			dn0 = min((max((dn0, bwmin0)), bwmax0))
			dn1 = min((max((dn1, bwmin1)), bwmax1))
		}
		n_cov = 0
		if (covadj) {
			if (rows(fuzzy) == rows(y) & bwparam == "main") outcomes = y, fuzzy
			else outcomes = y
			design = rd2d_basis1(u, p)
			w0full = rd2d_distance_weights(sdist, dn0, 0, kernel)
			w1full = rd2d_distance_weights(sdist, dn1, 1, kernel)
			gamma = rd2d_covariate_gamma(design, w0full, w1full, outcomes, covs)
			adjusted = rd2d_apply_covariates(outcomes, covs, gamma)
			ybase = adjusted[,1]
			if (rows(fuzzy) == rows(y) & bwparam == "main") fuzzybase = adjusted[,2]
			n_cov = rd2d_covariate_rank(design, w0full, w1full, covs)
		}
		y0 = select(ybase, idx0)
		y1 = select(ybase, idx1)
		thr0 = rd2d_quantile(u0, cqt)
		thr1 = rd2d_quantile(u1, cqt)
		fit0 = rd2d_poly_fit1(select(y0, u0 :<= thr0), select(u0, u0 :<= thr0), p + 1)
		fit1 = rd2d_poly_fit1(select(y1, u1 :<= thr1), select(u1, u1 :<= thr1), p + 1)

		if (rows(fuzzy) == rows(y) & bwparam == "main") {
			f0 = select(fuzzybase, idx0)
			f1 = select(fuzzybase, idx1)
			mu0 = rd2d_distance_intercepts_multi(y0, f0, u0, dn0, p, kernel)
			mu1 = rd2d_distance_intercepts_multi(y1, f1, u1, dn1, p, kernel)
			tauitt = mu1[1] - mu0[1]
			taufs = mu1[2] - mu0[2]
			if (taufs < . & abs(taufs) > sqrt(epsilon(1))) {
				graditt = 1 / taufs
				gradfs = -tauitt / taufs^2
				y0 = graditt :* y0 + gradfs :* f0
				y1 = graditt :* y1 + gradfs :* f1
				fit0 = rd2d_poly_fit1(select(y0, u0 :<= thr0), select(u0, u0 :<= thr0), p + 1)
				fit1 = rd2d_poly_fit1(select(y1, u1 :<= thr1), select(u1, u1 :<= thr1), p + 1)
			}
		}

		eN0 = sum(rd2d_kernel(u0:/dn0, kernel) :> 0)
		eN1 = sum(rd2d_kernel(u1:/dn1, kernel) :> 0)
		scales = rd2d_bw_vce_scales(fitmethod, vce, covadj, eN0, eN1, p + 1, n_cov)
		c0 = rd2d_distance_side_constants_adj(u0, y0, dn0, p, fit0, vce, kernel, scales[1])
		c1 = rd2d_distance_side_constants_adj(u1, y1, dn1, p, fit1, vce, kernel, scales[2])
		results[j,.] = (., ., c0[1], c1[1], c0[2], c1[2], c0[3], c1[3],
			c0[4], c1[4], bwmin0, bwmin1, bwmax0, bwmax1)
	}

	if (rd2d_bwbase(bwselect) == "mserd") {
		for (j=1; j<=neval; j++) {
			hn = (2 * (results[j,5] + results[j,6]) /
				((2*p + 2) * ((results[j,3] - results[j,4])^2 +
				scaleregul*results[j,7] + scaleregul*results[j,8])))^(1/(2*p + 4))
			if (bwcheck > 0) hn = min((max((hn, results[j,11], results[j,12])),
				max((results[j,13], results[j,14]))))
			results[j,1] = hn
			results[j,2] = hn
		}
	}
	else if (rd2d_bwbase(bwselect) == "msetwo") {
		for (j=1; j<=neval; j++) {
			h0 = (2 * results[j,5] / ((2*p + 2) *
				(results[j,3]^2 + scaleregul*results[j,7])))^(1/(2*p + 4))
			h1 = (2 * results[j,6] / ((2*p + 2) *
				(results[j,4]^2 + scaleregul*results[j,8])))^(1/(2*p + 4))
			if (bwcheck > 0) {
				h0 = min((max((h0, results[j,11])), results[j,13]))
				h1 = min((max((h1, results[j,12])), results[j,14]))
			}
			results[j,1] = h0
			results[j,2] = h1
		}
	}
	else if (rd2d_bwbase(bwselect) == "imserd") {
		VV = mean(results[,5]) + mean(results[,6])
		BB = mean((results[,3] :- results[,4]):^2 :+
			scaleregul:*results[,7] :+ scaleregul:*results[,8])
		hn = (2 * VV / ((2*p + 2) * BB))^(1/(2*p + 4))
		for (j=1; j<=neval; j++) {
			results[j,1] = hn
			results[j,2] = hn
			if (bwcheck > 0) {
				results[j,1] = min((max((results[j,1], results[j,11], results[j,12])),
					max((results[j,13], results[j,14]))))
				results[j,2] = results[j,1]
			}
		}
	}
	else {
		BB0 = mean(results[,3]:^2 :+ scaleregul:*results[,7])
		BB1 = mean(results[,4]:^2 :+ scaleregul:*results[,8])
		h0 = (2 * mean(results[,5]) / ((2*p + 2) * BB0))^(1/(2*p + 4))
		h1 = (2 * mean(results[,6]) / ((2*p + 2) * BB1))^(1/(2*p + 4))
		for (j=1; j<=neval; j++) {
			results[j,1] = h0
			results[j,2] = h1
			if (bwcheck > 0) {
				results[j,1] = min((max((results[j,1], results[j,11])), results[j,13]))
				results[j,2] = min((max((results[j,2], results[j,12])), results[j,14]))
			}
		}
	}

	smooth = rd2d_is_cer(bwselect) ? 1/(p+4) : 1/(2*p+4)
	if (rd2d_is_cer(bwselect)) {
		for (j=1; j<=neval; j++) {
			sdist = dist[,j]
			M0 = rows(uniqrows(select(sdist, sdist:<0)))
			M1 = rows(uniqrows(select(sdist, sdist:>=0)))
			M = M0 + M1
			if (rd2d_is_common(bwselect)) {
				factor = rd2d_cer_factor(M, p)
				results[j,1] = results[j,1] * factor
				results[j,2] = results[j,2] * factor
			}
			else {
				results[j,1] = results[j,1] * rd2d_cer_factor(M0, p)
				results[j,2] = results[j,2] * rd2d_cer_factor(M1, p)
			}
		}
	}
	if (sum(kinkpos) > 0) {
		distkink = rd2d_distance_to_known_kink(b, kinkpos)
		for (j=1; j<=neval; j++) {
			sdist = dist[,j]
			M0 = max((sum(sdist:<0), 1))
			M1 = max((sum(sdist:>=0), 1))
			M = M0 + M1
			if (rd2d_is_common(bwselect)) {
				factor = M^(smooth - 1/4)
				h0 = results[j,1] * factor
				h0 = max((h0, distkink[j]))
				results[j,1] = min((results[j,1], h0))
				results[j,2] = min((results[j,2], h0))
			}
			else {
				h0 = results[j,1] * M0^(smooth - 1/4)
				h1 = results[j,2] * M1^(smooth - 1/4)
				results[j,1] = min((results[j,1], max((h0, distkink[j]))))
				results[j,2] = min((results[j,2], max((h1, distkink[j]))))
			}
		}
	}
	if (kink1) {
		for (j=1; j<=neval; j++) {
			sdist = dist[,j]
			M0 = max((sum(sdist:<0), 1))
			M1 = max((sum(sdist:>=0), 1))
			M = M0 + M1
			if (rd2d_is_common(bwselect)) {
				factor = M^(smooth - 1/4)
				results[j,1] = results[j,1] * factor
				results[j,2] = results[j,2] * factor
			}
			else {
				results[j,1] = results[j,1] * M0^(smooth - 1/4)
				results[j,2] = results[j,2] * M1^(smooth - 1/4)
			}
		}
	}

	H = results[,1..2]
	bw = b, H, results[,9], results[,10]
	return(bw)
}

real matrix rd2d_apply_distance_bwcheck(real matrix H, real matrix dist,
	real scalar bwcheck)
{
	real matrix out
	real colvector sdist, u0, u1
	real rowvector b0, b1
	real scalar j

	out = H
	if (bwcheck <= 0) return(out)
	for (j=1; j<=cols(dist); j++) {
		sdist = dist[,j]
		u0 = abs(select(sdist, sdist:<0))
		u1 = abs(select(sdist, sdist:>=0))
		b0 = rd2d_kth_bounds(u0, bwcheck)
		b1 = rd2d_kth_bounds(u1, bwcheck)
		out[j,1] = min((max((out[j,1], b0[1])), b0[2]))
		out[j,2] = min((max((out[j,2], b1[1])), b1[2]))
	}
	return(out)
}

struct rd2d_orderfit scalar rd2d_fit_distance_order(real matrix outcomes, real matrix dist,
	real matrix H, real scalar p, string scalar kernel, string scalar vce,
	real colvector cluster, string scalar fitmethod, real matrix covs)
{
	struct rd2d_orderfit scalar out
	struct rd2d_lfit scalar fit0, fit1
	real matrix design, design0, design1, outcomes_fit, outcomes0, outcomes1, gamma
	real colvector sdist, u, u0, u1, w0, w1, active, clactive, cl0, cl1, cluster_groups
	real colvector idx0, idx1
	real colvector side0, side1, cluster0, cluster1
	real rowvector target
	real scalar j, n, neval, nout, ninf, covadj, n_cov, clustered, sidefit, e0, e1, kpoly
	real scalar scale0, scale1
	string scalar vce_local

	n = rows(dist)
	neval = cols(dist)
	nout = cols(outcomes)
	covadj = (cols(covs) > 0)
	clustered = (rows(cluster) == n)
	sidefit = (!covadj & (clustered | neval <= 6))
	cluster_groups = clustered ? uniqrows(sort(cluster,1)) : J(0,1,.)
	ninf = (clustered ? rows(cluster_groups) : n)
	fitmethod = strlower(strtrim(fitmethod))
	vce = strlower(strtrim(vce))
	out.mu0 = J(neval, nout, .)
	out.mu1 = J(neval, nout, .)
	out.se0 = J(neval, nout, .)
	out.se1 = J(neval, nout, .)
	out.N0 = J(neval, 1, 0)
	out.N1 = J(neval, 1, 0)
	out.infl0_1 = J(neval, ninf, 0)
	out.infl1_1 = J(neval, ninf, 0)
	out.infl0_2 = J(neval, ninf, 0)
	out.infl1_2 = J(neval, ninf, 0)
	target = rd2d_target1(p)

	for (j=1; j<=neval; j++) {
		sdist = dist[,j]
		n_cov = 0
		if (sidefit) {
			side0 = (sdist :< 0)
			side1 = (sdist :>= 0)
			idx0 = selectindex(side0)
			idx1 = selectindex(side1)
			u0 = abs(sdist[idx0])
			u1 = abs(sdist[idx1])
			design0 = rd2d_basis1(u0, p)
			design1 = rd2d_basis1(u1, p)
			w0 = rd2d_kernel_nonnegative(u0:/H[j,1], kernel) :/ max((H[j,1]^2, epsilon(1)))
			w1 = rd2d_kernel_nonnegative(u1:/H[j,2], kernel) :/ max((H[j,2]^2, epsilon(1)))
			outcomes0 = outcomes[idx0,.]
			outcomes1 = outcomes[idx1,.]
			cluster0 = clustered ? cluster[idx0] : J(0,1,.)
			cluster1 = clustered ? cluster[idx1] : J(0,1,.)
			e0 = sum(w0 :> 0)
			e1 = sum(w1 :> 0)
			kpoly = cols(design0)
			vce_local = ((fitmethod == "joint" | covadj) & vce == "hc1") ? "hc0" : vce
			scale0 = .
			scale1 = .
			if ((fitmethod == "joint" | covadj) & (vce == "hc1" | clustered)) {
				if (fitmethod == "joint") {
					clactive = select(cluster0, w0 :> 0) \ select(cluster1, w1 :> 0)
					scale0 = rd2d_joint_vce_scale(vce, e0 + e1, 2*kpoly + n_cov, clactive, clustered)
					scale1 = scale0
				}
				else {
					cl0 = select(cluster0, w0 :> 0)
					cl1 = select(cluster1, w1 :> 0)
					scale0 = rd2d_joint_vce_scale(vce, e0, kpoly + n_cov, cl0, clustered)
					scale1 = rd2d_joint_vce_scale(vce, e1, kpoly + n_cov, cl1, clustered)
				}
			}
			fit0 = rd2d_local_fit_scaled(design0, w0, outcomes0, target, vce_local, cluster0, scale0, cluster_groups)
			fit1 = rd2d_local_fit_scaled(design1, w1, outcomes1, target, vce_local, cluster1, scale1, cluster_groups)
		}
		else {
			u = abs(sdist)
			design = rd2d_basis1(u, p)
			w0 = rd2d_distance_weights(sdist, H[j,1], 0, kernel)
			w1 = rd2d_distance_weights(sdist, H[j,2], 1, kernel)
			outcomes_fit = outcomes
			if (covadj) {
			gamma = rd2d_covariate_gamma(design, w0, w1, outcomes, covs)
			n_cov = rd2d_covariate_rank(design, w0, w1, covs)
			outcomes_fit = rd2d_apply_covariates(outcomes, covs, gamma)
			}
			e0 = sum(w0 :> 0)
			e1 = sum(w1 :> 0)
			kpoly = cols(design)
			vce_local = ((fitmethod == "joint" | covadj) & vce == "hc1") ? "hc0" : vce
			scale0 = .
			scale1 = .
			if ((fitmethod == "joint" | covadj) & (vce == "hc1" | clustered)) {
				if (fitmethod == "joint") {
					active = (w0 :> 0) :| (w1 :> 0)
					clactive = clustered ? select(cluster, active) : J(0,1,.)
					scale0 = rd2d_joint_vce_scale(vce, e0 + e1, 2*kpoly + n_cov, clactive, clustered)
					scale1 = scale0
				}
				else {
					cl0 = clustered ? select(cluster, w0 :> 0) : J(0,1,.)
					cl1 = clustered ? select(cluster, w1 :> 0) : J(0,1,.)
					scale0 = rd2d_joint_vce_scale(vce, e0, kpoly + n_cov, cl0, clustered)
					scale1 = rd2d_joint_vce_scale(vce, e1, kpoly + n_cov, cl1, clustered)
				}
			}
			fit0 = rd2d_local_fit_scaled(design, w0, outcomes_fit, target, vce_local, cluster, scale0, cluster_groups)
			fit1 = rd2d_local_fit_scaled(design, w1, outcomes_fit, target, vce_local, cluster, scale1, cluster_groups)
		}
		out.mu0[j,.] = fit0.estimate
		out.mu1[j,.] = fit1.estimate
		out.se0[j,.] = fit0.se
		out.se1[j,.] = fit1.se
		out.N0[j] = fit0.n_eff
		out.N1[j] = fit1.n_eff
		if (sidefit & !clustered) {
			out.infl0_1[j,idx0] = fit0.influence[,1]'
			out.infl1_1[j,idx1] = fit1.influence[,1]'
			if (nout >= 2) {
				out.infl0_2[j,idx0] = fit0.influence[,2]'
				out.infl1_2[j,idx1] = fit1.influence[,2]'
			}
		}
		else {
			out.infl0_1[j,.] = fit0.influence[,1]'
			out.infl1_1[j,.] = fit1.influence[,1]'
			if (nout >= 2) {
				out.infl0_2[j,.] = fit0.influence[,2]'
				out.infl1_2[j,.] = fit1.influence[,2]'
			}
		}
	}
	return(out)
}

void rd2d_distance_stata(string scalar main_name, string scalar bw_name,
	string scalar main0_name, string scalar main1_name, string scalar itt_name,
	string scalar fs_name, string scalar itt0_name, string scalar itt1_name,
	string scalar fs0_name, string scalar fs1_name, string scalar vmain_name,
	string scalar vitt_name, string scalar vfs_name, string scalar vitt0_name,
	string scalar vitt1_name, string scalar vfs0_name, string scalar vfs1_name,
	string scalar yvar,
	string scalar distvars, string scalar fuzzyvar, string scalar clustvar,
	string scalar touse, string scalar bstr, string scalar hstr, real scalar p,
	real scalar q, string scalar kernel, string scalar vce, real scalar level,
	string scalar side, string scalar bwselect, string scalar bwparam,
	real scalar bwcheck, real scalar kink1, real scalar kink2,
	string scalar kinkposstr, real scalar scaleregul, real scalar cqt,
	string scalar fitmethod, string scalar covsvars)
{
	real matrix data, dist, b, bw, H, Hr, outcomes, covs
	real matrix main, main0, main1, itt, fs, itt0, itt1, fs0, fs1
	real matrix Vmain, Vitt, Vfs, Vitt0, Vitt1, Vfs0, Vfs1
	real colvector y, fuzzy, cluster, taup, tauq, sep, seq
	real colvector itt_p, itt_q, fs_p, fs_q, se_itt_p, se_itt_q, se_fs_p, se_fs_q
	real matrix inflp, inflq, infl_itt_p, infl_itt_q, infl_fs_p, infl_fs_q
	struct rd2d_orderfit scalar fitp, fitq
	real scalar neval, smooth, covstart, covend

	data = st_data(., yvar+" "+distvars+
		(fuzzyvar=="" ? "" : " "+fuzzyvar)+
		(covsvars=="" ? "" : " "+covsvars)+
		(clustvar=="" ? "" : " "+clustvar), touse)
	y = data[,1]
	neval = cols(tokens(distvars))
	dist = data[,2..(1+neval)]
	if (bstr == "") b = J(neval, 2, .)
	else b = rd2d_parse_numlist(bstr, 2)
	if (rows(b) != neval) {
		errprintf("b() must contain two values per distance variable\n")
		_error(198)
	}
	if (fuzzyvar != "") fuzzy = data[,2+neval]
	else fuzzy = J(0,1,.)
	if (covsvars != "") {
		covstart = 2 + neval + (fuzzyvar != "")
		covend = covstart + cols(tokens(covsvars)) - 1
		covs = data[,covstart..covend]
	}
	else covs = J(rows(data), 0, .)
	if (clustvar != "") cluster = data[,cols(data)]
	else cluster = J(0,1,.)
	fitmethod = (strlower(strtrim(fitmethod)) == "separate" ? "separate" : "joint")
	if (q < 0) q = kink1 ? p : p + 1

	bw = rd2d_distance_bw_full(dist, b, hstr, y, fuzzy, p, kernel,
		bwselect, bwparam, vce, bwcheck, kink1, kinkposstr, scaleregul, cqt,
		fitmethod, covs)
	H = bw[,(3,4)]
	H = rd2d_apply_distance_bwcheck(H, dist, bwcheck)
	Hr = H
	if (kink2) {
		smooth = rd2d_is_cer(bwselect) ? 1/(p+4) : 1/(2*p+4)
		if (rd2d_is_common(bwselect)) Hr = H :* (rows(dist)^(smooth-1/4))
		else {
			Hr = H
			Hr[,1] = H[,1] :* (rows(dist)^(smooth-1/4))
			Hr[,2] = H[,2] :* (rows(dist)^(smooth-1/4))
		}
		Hr = rd2d_apply_distance_bwcheck(Hr, dist, bwcheck)
	}
	if (fuzzyvar == "") outcomes = y
	else outcomes = y, fuzzy
	fitp = rd2d_fit_distance_order(outcomes, dist, H, p, kernel, vce, cluster, fitmethod, covs)
	fitq = (q == p & all(vec(H :== Hr)) ? fitp : rd2d_fit_distance_order(outcomes, dist, Hr, q, kernel, vce, cluster, fitmethod, covs))

	if (fuzzyvar == "") {
		taup = fitp.mu1[,1] - fitp.mu0[,1]
		tauq = fitq.mu1[,1] - fitq.mu0[,1]
		if (rows(cluster) == rows(y) & fitmethod == "separate") {
			inflp = fitp.infl1_1, -fitp.infl0_1
			inflq = fitq.infl1_1, -fitq.infl0_1
		}
		else {
			inflp = fitp.infl1_1 - fitp.infl0_1
			inflq = fitq.infl1_1 - fitq.infl0_1
		}
		sep = sqrt(diagonal(inflp*inflp'))
		seq = sqrt(diagonal(inflq*inflq'))
		main = rd2d_ci_table(b, taup, sep, tauq, seq, (H,Hr), fitp.N0, fitp.N1, level, side, 1)
		main0 = rd2d_ci_table(b, fitp.mu0[,1], fitp.se0[,1], fitq.mu0[,1], fitq.se0[,1],
			(H[,1], J(rows(H),1,.), Hr[,1], J(rows(H),1,.)), fitp.N0, J(rows(H),1,.), level, side, 1)
		main1 = rd2d_ci_table(b, fitp.mu1[,1], fitp.se1[,1], fitq.mu1[,1], fitq.se1[,1],
			(J(rows(H),1,.), H[,2], J(rows(H),1,.), Hr[,2]), J(rows(H),1,.), fitp.N1, level, side, 1)
		itt = J(0,0,.)
		fs = J(0,0,.)
		itt0 = J(0,0,.)
		itt1 = J(0,0,.)
		fs0 = J(0,0,.)
		fs1 = J(0,0,.)
		Vmain = inflq * inflq'
		Vitt = J(0,0,.)
		Vfs = J(0,0,.)
		Vitt0 = J(0,0,.)
		Vitt1 = J(0,0,.)
		Vfs0 = J(0,0,.)
		Vfs1 = J(0,0,.)
	}
	else {
		itt_p = fitp.mu1[,1] - fitp.mu0[,1]
		itt_q = fitq.mu1[,1] - fitq.mu0[,1]
		fs_p = fitp.mu1[,2] - fitp.mu0[,2]
		fs_q = fitq.mu1[,2] - fitq.mu0[,2]
		taup = itt_p :/ fs_p
		tauq = itt_q :/ fs_q
		if (rows(cluster) == rows(y) & fitmethod == "separate") {
			infl_itt_p = fitp.infl1_1, -fitp.infl0_1
			infl_itt_q = fitq.infl1_1, -fitq.infl0_1
			infl_fs_p = fitp.infl1_2, -fitp.infl0_2
			infl_fs_q = fitq.infl1_2, -fitq.infl0_2
		}
		else {
			infl_itt_p = fitp.infl1_1 - fitp.infl0_1
			infl_itt_q = fitq.infl1_1 - fitq.infl0_1
			infl_fs_p = fitp.infl1_2 - fitp.infl0_2
			infl_fs_q = fitq.infl1_2 - fitq.infl0_2
		}
		inflp = (infl_itt_p :/ (fs_p*J(1,cols(infl_itt_p),1))) :-
			((itt_p:/fs_p:^2)*J(1,cols(infl_fs_p),1)) :* infl_fs_p
		inflq = (infl_itt_q :/ (fs_q*J(1,cols(infl_itt_q),1))) :-
			((itt_q:/fs_q:^2)*J(1,cols(infl_fs_q),1)) :* infl_fs_q
		sep = sqrt(diagonal(inflp*inflp'))
		seq = sqrt(diagonal(inflq*inflq'))
		se_itt_p = sqrt(diagonal(infl_itt_p*infl_itt_p'))
		se_itt_q = sqrt(diagonal(infl_itt_q*infl_itt_q'))
		se_fs_p = sqrt(diagonal(infl_fs_p*infl_fs_p'))
		se_fs_q = sqrt(diagonal(infl_fs_q*infl_fs_q'))
		main = rd2d_ci_table(b, taup, sep, tauq, seq, (H,Hr), fitp.N0, fitp.N1, level, side, 1)
		itt = rd2d_ci_table(b, itt_p, se_itt_p, itt_q, se_itt_q, (H,Hr), fitp.N0, fitp.N1, level, side, 1)
		fs = rd2d_ci_table(b, fs_p, se_fs_p, fs_q, se_fs_q, (H,Hr), fitp.N0, fitp.N1, level, side, 1)
		itt0 = rd2d_ci_table(b, fitp.mu0[,1], fitp.se0[,1], fitq.mu0[,1], fitq.se0[,1],
			(H[,1], J(rows(H),1,.), Hr[,1], J(rows(H),1,.)), fitp.N0, J(rows(H),1,.), level, side, 1)
		itt1 = rd2d_ci_table(b, fitp.mu1[,1], fitp.se1[,1], fitq.mu1[,1], fitq.se1[,1],
			(J(rows(H),1,.), H[,2], J(rows(H),1,.), Hr[,2]), J(rows(H),1,.), fitp.N1, level, side, 1)
		fs0 = rd2d_ci_table(b, fitp.mu0[,2], fitp.se0[,2], fitq.mu0[,2], fitq.se0[,2],
			(H[,1], J(rows(H),1,.), Hr[,1], J(rows(H),1,.)), fitp.N0, J(rows(H),1,.), level, side, 1)
		fs1 = rd2d_ci_table(b, fitp.mu1[,2], fitp.se1[,2], fitq.mu1[,2], fitq.se1[,2],
			(J(rows(H),1,.), H[,2], J(rows(H),1,.), Hr[,2]), J(rows(H),1,.), fitp.N1, level, side, 1)
		main0 = J(0,0,.)
		main1 = J(0,0,.)
		Vmain = inflq * inflq'
		Vitt = infl_itt_q * infl_itt_q'
		Vfs = infl_fs_q * infl_fs_q'
		Vitt0 = fitq.infl0_1 * fitq.infl0_1'
		Vitt1 = fitq.infl1_1 * fitq.infl1_1'
		Vfs0 = fitq.infl0_2 * fitq.infl0_2'
		Vfs1 = fitq.infl1_2 * fitq.infl1_2'
	}
	bw = b, H, fitp.N0, fitp.N1
	st_matrix(main_name, main)
	st_matrix(bw_name, bw)
	st_matrix(main0_name, main0)
	st_matrix(main1_name, main1)
	st_matrix(itt_name, itt)
	st_matrix(fs_name, fs)
	st_matrix(itt0_name, itt0)
	st_matrix(itt1_name, itt1)
	st_matrix(fs0_name, fs0)
	st_matrix(fs1_name, fs1)
	st_matrix(vmain_name, Vmain)
	st_matrix(vitt_name, Vitt)
	st_matrix(vfs_name, Vfs)
	st_matrix(vitt0_name, Vitt0)
	st_matrix(vitt1_name, Vitt1)
	st_matrix(vfs0_name, Vfs0)
	st_matrix(vfs1_name, Vfs1)
}

void rdbw2d_distance_stata(string scalar bw_name, string scalar yvar,
	string scalar distvars, string scalar fuzzyvar, string scalar clustvar,
	string scalar touse, string scalar bstr, real scalar p, string scalar kernel,
	string scalar bwselect, string scalar bwparam, string scalar vce,
	real scalar bwcheck, real scalar kink1, string scalar kinkposstr,
	real scalar scaleregul, real scalar cqt, string scalar fitmethod,
	string scalar covsvars)
{
	real matrix data, dist, b, bw, covs
	real colvector y, fuzzy
	real scalar neval, covstart, covend

	data = st_data(., yvar+" "+distvars+
		(fuzzyvar=="" ? "" : " "+fuzzyvar)+
		(covsvars=="" ? "" : " "+covsvars)+
		(clustvar=="" ? "" : " "+clustvar), touse)
	y = data[,1]
	neval = cols(tokens(distvars))
	dist = data[,2..(1+neval)]
	if (bstr == "") b = J(neval, 2, .)
	else b = rd2d_parse_numlist(bstr, 2)
	if (rows(b) != neval) {
		errprintf("b() must contain two values per distance variable\n")
		_error(198)
	}
	if (fuzzyvar != "") fuzzy = data[,2+neval]
	else fuzzy = J(0,1,.)
	if (covsvars != "") {
		covstart = 2 + neval + (fuzzyvar != "")
		covend = covstart + cols(tokens(covsvars)) - 1
		covs = data[,covstart..covend]
	}
	else covs = J(rows(data), 0, .)
	fitmethod = (strlower(strtrim(fitmethod)) == "separate" ? "separate" : "joint")
	bw = rd2d_distance_bw_full(dist, b, "", y, fuzzy, p, kernel,
		bwselect, bwparam, vce, bwcheck, kink1, kinkposstr, scaleregul, cqt,
		fitmethod, covs)
	st_matrix(bw_name, bw)
}

end

if "`rd2d_loadonly'" != "1" {
	capture erase lrd2d.mlib
	mata: mata mlib create lrd2d, replace
	mata: mata mlib add lrd2d rd2d_lfit() rd2d_orderfit() ///
		rd2d_bwfit() rd2d_bwconst() rd2d_polyfit() rd2d_covcomp()
	mata: mata mlib add lrd2d rd2d_mlib_loaded() rd2d_cluster_sums() ///
		rd2d_parse_numlist() ///
		rd2d_parse_kink_position() rd2d_distance_to_known_kink() ///
		rd2d_fact() rd2d_cer_factor() rd2d_bwbase() rd2d_is_cer() ///
		rd2d_is_common() rd2d_kernel() rd2d_kernel_nonnegative() ///
		rd2d_basis1() rd2d_basis2() ///
		rd2d_target2() rd2d_target1() rd2d_location_weights() ///
		rd2d_distance_weights() rd2d_local_fit() rd2d_local_fit_scaled() ///
		rd2d_joint_vce_scale() rd2d_projection_components() ///
		rd2d_covariate_gamma() rd2d_covariate_rank() ///
		rd2d_apply_covariates() rd2d_bw_vce_scales()
	mata: mata mlib add lrd2d rd2d_basis_count2() rd2d_H2() ///
		rd2d_invH2() rd2d_vce_const2() rd2d_lm_exact2() ///
		rd2d_lm_multi_exact2() rd2d_get_coeff_exact2() ///
		rd2d_bw_v2_exact2_adj() rd2d_bw_v2_exact2() ///
		rd2d_rot_location() rd2d_rot_distance() ///
		rd2d_kth_bounds() rd2d_quantile() rd2d_apply_location_bwcheck() ///
		rd2d_location_bw_simple() rd2d_location_bw_full()
	mata: mata mlib add lrd2d rd2d_fit_location_order() rd2d_ci_table() ///
		rd2d_location_stata() rdbw2d_location_stata() ///
		rd2d_distance_bw_simple() rd2d_poly_fit1() rd2d_vce_const1() ///
		rd2d_distance_intercepts_multi() rd2d_distance_side_constants() ///
		rd2d_distance_side_constants_adj() rd2d_distance_bw_full() ///
		rd2d_apply_distance_bwcheck() ///
		rd2d_fit_distance_order() rd2d_distance_stata() ///
		rdbw2d_distance_stata()
	mata: mata mlib index
}
