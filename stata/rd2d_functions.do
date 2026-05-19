********************************************************************************
* RD2D STATA PACKAGE -- Mata functions
* Authors: Matias D. Cattaneo, Rocio Titiunik, Ruiqi Rae Yu
********************************************************************************
*!version 0.1.0  2026-05-19

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
	struct rd2d_lfit scalar out
	real colvector idx, w, leverage, adj, score, sqrtadj, ckeep, groups
	real matrix X, Y, gram, invG, beta, fitted, resid, Xw, infl_kept, infl_full, cov, part
	real rowvector row
	real scalar n, k, m, i, g, scale, has_cluster

	n = rows(design)
	k = cols(design)
	m = cols(outcomes)
	has_cluster = (rows(cluster) == n)
	out.estimate = J(1, m, .)
	out.se = J(1, m, .)
	out.n_eff = 0
	out.influence = J(has_cluster ? rows(uniqrows(sort(cluster,1))) : n, m, 0)

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

	leverage = rowsum((X * invG) :* X) :* w
	adj = J(n,1,1)
	vce = strlower(strtrim(vce))
	if (vce == "hc2") adj = 1 :/ rowmax((1 :- leverage, J(n,1,1e-8)))
	if (vce == "hc3") adj = (1 :/ rowmax((1 :- leverage, J(n,1,1e-8)))):^2

	row = target * invG
	score = Xw * row'
	sqrtadj = sqrt(adj)
	infl_kept = (score * J(1,m,1)) :* resid :* (sqrtadj * J(1,m,1))

	scale = 1
	if (has_cluster) {
		groups = uniqrows(sort(cluster,1))
		ckeep = cluster[idx]
		infl_full = J(rows(groups), m, 0)
		for (g=1; g<=rows(groups); g++) {
			part = select(infl_kept, ckeep :== groups[g])
			if (rows(part) > 0) infl_full[g,.] = colsum(part)
		}
		if (vce == "hc1" & rows(groups) > 1 & n > k) {
			scale = (rows(groups)/(rows(groups)-1)) * ((n-1)/(n-k))
		}
		cov = scale * infl_full' * infl_full
		out.influence = sqrt(scale) * infl_full
	}
	else {
		if (vce == "hc1" & n > k) scale = n/(n-k)
		cov = scale * infl_kept' * infl_kept
		infl_full = J(rows(design), m, 0)
		infl_full[idx,.] = sqrt(scale) * infl_kept
		out.influence = infl_full
	}
	out.estimate = target * beta
	out.se = sqrt(diagonal(cov))'
	return(out)
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

struct rd2d_orderfit scalar rd2d_fit_location_order(real matrix x, real colvector d,
	real matrix b, real matrix outcomes, real matrix H, real scalar p,
	real rowvector deriv, real matrix tangvec, string scalar kernel,
	string scalar kerneltype, string scalar vce, real colvector cluster)
{
	struct rd2d_orderfit scalar out
	struct rd2d_lfit scalar fit0, fit1
	real matrix centered, design
	real colvector w0, w1
	real rowvector target
	real scalar j, n, neval, nout, ninf

	n = rows(x)
	neval = rows(b)
	nout = cols(outcomes)
	ninf = (rows(cluster)==n ? rows(uniqrows(sort(cluster,1))) : n)
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

	for (j=1; j<=neval; j++) {
		centered = x :- (J(n,1,1)*b[j,.])
		design = rd2d_basis2(centered, p)
		target = rd2d_target2(p, deriv, tangvec, j)
		w0 = rd2d_location_weights(centered, H[j,(1,2)], kernel, kerneltype) :* (d:==0)
		w1 = rd2d_location_weights(centered, H[j,(3,4)], kernel, kerneltype) :* (d:!=0)
		fit0 = rd2d_local_fit(design, w0, outcomes, target, vce, cluster)
		fit1 = rd2d_local_fit(design, w1, outcomes, target, vce, cluster)
		out.mu0[j,.] = fit0.estimate
		out.mu1[j,.] = fit1.estimate
		out.se0[j,.] = fit0.se
		out.se1[j,.] = fit1.se
		out.N0[j] = fit0.n_eff
		out.N1[j] = fit1.n_eff
		out.infl0_1[j,.] = fit0.influence[,1]'
		out.infl1_1[j,.] = fit1.influence[,1]'
		if (nout >= 2) {
			out.infl0_2[j,.] = fit0.influence[,2]'
			out.infl1_2[j,.] = fit1.influence[,2]'
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
	string scalar side, string scalar bwselect, real scalar bwcheck,
	real scalar stdvars)
{
	real matrix data, x, b, bw, H, outcomes
	real matrix main, main0, main1, itt, fs, itt0, itt1, fs0, fs1
	real matrix Vmain, Vitt, Vfs, Vitt0, Vitt1, Vfs0, Vfs1
	real colvector y, d, fuzzy, cluster, estp, estq, sep, seq, taup, tauq
	real colvector itt_p, itt_q, fs_p, fs_q, se_itt_p, se_itt_q, se_fs_p, se_fs_q
	real matrix inflp, inflq, infl_itt_p, infl_itt_q, infl_fs_p, infl_fs_q
	real matrix deriv_raw
	real rowvector deriv
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
		(clustvar=="" ? "" : " "+clustvar), touse)
	y = data[,1]
	x = data[,(2,3)]
	d = data[,4]
	if (fuzzyvar != "") fuzzy = data[,5]
	else fuzzy = J(0,1,.)
	if (clustvar != "") cluster = data[,cols(data)]
	else cluster = J(0,1,.)

	bw = rd2d_location_bw_simple(y, x, d, b, hstr, p, kernel, kerneltype,
		bwselect, bwcheck, stdvars)
	H = bw[,(3,4,5,6)]
	if (fuzzyvar == "") outcomes = y
	else outcomes = y, fuzzy
	fitp = rd2d_fit_location_order(x, d, b, outcomes, H, p, deriv, tangvec, kernel, kerneltype, vce, cluster)
	fitq = (q == p ? fitp : rd2d_fit_location_order(x, d, b, outcomes, H, q, deriv, tangvec, kernel, kerneltype, vce, cluster))

	if (fuzzyvar == "") {
		taup = fitp.mu1[,1] - fitp.mu0[,1]
		tauq = fitq.mu1[,1] - fitq.mu0[,1]
		inflp = fitp.infl1_1 - fitp.infl0_1
		inflq = fitq.infl1_1 - fitq.infl0_1
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
		infl_itt_p = fitp.infl1_1 - fitp.infl0_1
		infl_itt_q = fitq.infl1_1 - fitq.infl0_1
		infl_fs_p = fitp.infl1_2 - fitp.infl0_2
		infl_fs_q = fitq.infl1_2 - fitq.infl0_2
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
	string scalar bstr, real scalar p, string scalar kernel,
	string scalar kerneltype, string scalar bwselect, real scalar bwcheck,
	real scalar stdvars)
{
	real matrix data, x, b, bw
	real colvector y, d

	b = rd2d_parse_numlist(bstr, 2)
	data = st_data(., yvar+" "+x1var+" "+x2var+" "+dvar+
		(fuzzyvar=="" ? "" : " "+fuzzyvar)+
		(clustvar=="" ? "" : " "+clustvar), touse)
	y = data[,1]
	x = data[,(2,3)]
	d = data[,4]
	bw = rd2d_location_bw_simple(y, x, d, b, "", p, kernel, kerneltype,
		bwselect, bwcheck, stdvars)
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

struct rd2d_orderfit scalar rd2d_fit_distance_order(real matrix outcomes, real matrix dist,
	real matrix H, real scalar p, string scalar kernel, string scalar vce,
	real colvector cluster)
{
	struct rd2d_orderfit scalar out
	struct rd2d_lfit scalar fit0, fit1
	real matrix design
	real colvector sdist, u, w0, w1
	real rowvector target
	real scalar j, n, neval, nout, ninf

	n = rows(dist)
	neval = cols(dist)
	nout = cols(outcomes)
	ninf = (rows(cluster)==n ? rows(uniqrows(sort(cluster,1))) : n)
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
		u = abs(sdist)
		design = rd2d_basis1(u, p)
		w0 = rd2d_distance_weights(sdist, H[j,1], 0, kernel)
		w1 = rd2d_distance_weights(sdist, H[j,2], 1, kernel)
		fit0 = rd2d_local_fit(design, w0, outcomes, target, vce, cluster)
		fit1 = rd2d_local_fit(design, w1, outcomes, target, vce, cluster)
		out.mu0[j,.] = fit0.estimate
		out.mu1[j,.] = fit1.estimate
		out.se0[j,.] = fit0.se
		out.se1[j,.] = fit1.se
		out.N0[j] = fit0.n_eff
		out.N1[j] = fit1.n_eff
		out.infl0_1[j,.] = fit0.influence[,1]'
		out.infl1_1[j,.] = fit1.influence[,1]'
		if (nout >= 2) {
			out.infl0_2[j,.] = fit0.influence[,2]'
			out.infl1_2[j,.] = fit1.influence[,2]'
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
	string scalar side, string scalar bwselect, real scalar bwcheck,
	real scalar kink1, real scalar kink2)
{
	real matrix data, dist, b, bw, H, Hr, outcomes
	real matrix main, main0, main1, itt, fs, itt0, itt1, fs0, fs1
	real matrix Vmain, Vitt, Vfs, Vitt0, Vitt1, Vfs0, Vfs1
	real colvector y, fuzzy, cluster, taup, tauq, sep, seq
	real colvector itt_p, itt_q, fs_p, fs_q, se_itt_p, se_itt_q, se_fs_p, se_fs_q
	real matrix inflp, inflq, infl_itt_p, infl_itt_q, infl_fs_p, infl_fs_q
	struct rd2d_orderfit scalar fitp, fitq
	real scalar neval

	data = st_data(., yvar+" "+distvars+
		(fuzzyvar=="" ? "" : " "+fuzzyvar)+
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
	if (clustvar != "") cluster = data[,cols(data)]
	else cluster = J(0,1,.)
	if (q < 0) q = kink1 ? p : p + 1

	bw = rd2d_distance_bw_simple(dist, b, hstr, p, kernel, bwselect, bwcheck, kink1)
	H = bw[,(3,4)]
	Hr = H
	if (kink2) {
		if (rd2d_is_common(bwselect)) Hr = H :* (rows(dist)^(1/(2*p+4)-1/4))
		else {
			Hr = H
			Hr[,1] = H[,1] :* (rows(dist)^(1/(2*p+4)-1/4))
			Hr[,2] = H[,2] :* (rows(dist)^(1/(2*p+4)-1/4))
		}
	}
	if (fuzzyvar == "") outcomes = y
	else outcomes = y, fuzzy
	fitp = rd2d_fit_distance_order(outcomes, dist, H, p, kernel, vce, cluster)
	fitq = (q == p & all(vec(H :== Hr)) ? fitp : rd2d_fit_distance_order(outcomes, dist, Hr, q, kernel, vce, cluster))

	if (fuzzyvar == "") {
		taup = fitp.mu1[,1] - fitp.mu0[,1]
		tauq = fitq.mu1[,1] - fitq.mu0[,1]
		inflp = fitp.infl1_1 - fitp.infl0_1
		inflq = fitq.infl1_1 - fitq.infl0_1
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
		infl_itt_p = fitp.infl1_1 - fitp.infl0_1
		infl_itt_q = fitq.infl1_1 - fitq.infl0_1
		infl_fs_p = fitp.infl1_2 - fitp.infl0_2
		infl_fs_q = fitq.infl1_2 - fitq.infl0_2
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
	string scalar bwselect, real scalar bwcheck, real scalar kink1)
{
	real matrix data, dist, b, bw
	real scalar neval

	data = st_data(., yvar+" "+distvars+
		(fuzzyvar=="" ? "" : " "+fuzzyvar)+
		(clustvar=="" ? "" : " "+clustvar), touse)
	neval = cols(tokens(distvars))
	dist = data[,2..(1+neval)]
	if (bstr == "") b = J(neval, 2, .)
	else b = rd2d_parse_numlist(bstr, 2)
	if (rows(b) != neval) {
		errprintf("b() must contain two values per distance variable\n")
		_error(198)
	}
	bw = rd2d_distance_bw_simple(dist, b, "", p, kernel, bwselect, bwcheck, kink1)
	st_matrix(bw_name, bw)
}

end

if "`rd2d_loadonly'" != "1" {
	mata: mata mosave rd2d_parse_numlist(), replace
	mata: mata mosave rd2d_fact(), replace
	mata: mata mosave rd2d_cer_factor(), replace
	mata: mata mosave rd2d_bwbase(), replace
	mata: mata mosave rd2d_is_cer(), replace
	mata: mata mosave rd2d_is_common(), replace
	mata: mata mosave rd2d_kernel(), replace
	mata: mata mosave rd2d_basis1(), replace
	mata: mata mosave rd2d_basis2(), replace
	mata: mata mosave rd2d_target2(), replace
	mata: mata mosave rd2d_target1(), replace
	mata: mata mosave rd2d_location_weights(), replace
	mata: mata mosave rd2d_distance_weights(), replace
	mata: mata mosave rd2d_local_fit(), replace
	mata: mata mosave rd2d_rot_location(), replace
	mata: mata mosave rd2d_rot_distance(), replace
	mata: mata mosave rd2d_kth_bounds(), replace
	mata: mata mosave rd2d_apply_location_bwcheck(), replace
	mata: mata mosave rd2d_location_bw_simple(), replace
	mata: mata mosave rd2d_fit_location_order(), replace
	mata: mata mosave rd2d_ci_table(), replace
	mata: mata mosave rd2d_location_stata(), replace
	mata: mata mosave rdbw2d_location_stata(), replace
	mata: mata mosave rd2d_distance_bw_simple(), replace
	mata: mata mosave rd2d_fit_distance_order(), replace
	mata: mata mosave rd2d_distance_stata(), replace
	mata: mata mosave rdbw2d_distance_stata(), replace
}
