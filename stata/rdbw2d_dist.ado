********************************************************************************
* RD2D STATA PACKAGE -- rdbw2d_dist
* Authors: Matias D. Cattaneo, Rocio Titiunik, Ruiqi Rae Yu
********************************************************************************
*!version 0.2.0  2026-05-23

capture program drop rdbw2d_dist
program define rdbw2d_dist, eclass
	version 16.0
	syntax varlist(min=2 numeric) [if] [in] ///
		[, B(string asis) Fuzzy(varname numeric) P(integer 1) ///
		   KINKUNKnown(string) KINKPOSition(numlist) KERnel(string) ///
		   BWSELECT(string) BWPARAM(string) VCE(string) BWCHECK(integer -1) ///
		   MASSpoints(string) CLuster(varname numeric) SCALEREGUL(real 1) ///
		   CQT(real .5) ]

	marksample touse
	gettoken y distvars : varlist
	markout `touse' `y' `distvars'
	if "`fuzzy'" != "" markout `touse' `fuzzy'
	if "`cluster'" != "" markout `touse' `cluster'

	if "`kernel'" == "" local kernel "tri"
	if "`bwselect'" == "" local bwselect "mserd"
	if "`bwparam'" == "" local bwparam "main"
	if "`vce'" == "" local vce "hc1"
	if "`masspoints'" == "" local masspoints "check"
	if `bwcheck' < 0 local bwcheck = 20 + `p' + 1

	local kink1 = 0
	local ku = lower(strtrim("`kinkunknown'"))
	if inlist("`ku'", "on", "yes", "true", "1") local kink1 = 1
	else if "`ku'" != "" & !inlist("`ku'", "off", "no", "false", "0") {
		tokenize "`ku'"
		local kink1 = real("`1'")
	}
	if missing(`kink1') | !inlist(`kink1', 0, 1) {
		di as error "kinkunknown() must be on, off, or a 0/1 indicator"
		exit 198
	}
	if "`kinkposition'" != "" & `kink1' {
		di as error "use either kinkposition() or kinkunknown(), not both"
		exit 198
	}

	if "$RD2D_MATA_LOADED" != "1" {
		tempname rd2d_mlib_ok
		capture quietly mata: mata mlib index
		capture quietly mata: st_numscalar("`rd2d_mlib_ok'", rd2d_mlib_loaded())
		if _rc {
			local rd2d_loadonly 1
			quietly findfile rd2d_functions.do
			quietly do "`r(fn)'"
		}
		global RD2D_MATA_LOADED 1
	}

	tempname bw
	mata: rdbw2d_distance_stata("`bw'", "`y'", "`distvars'", "`fuzzy'", ///
		"`cluster'", "`touse'", st_local("b"), `p', "`kernel'", ///
		"`bwselect'", "`bwparam'", "`vce'", `bwcheck', `kink1', ///
		st_local("kinkposition"), `scaleregul', `cqt')

	matrix colnames `bw' = b1 b2 h0 h1 N_Co N_Tr
	ereturn clear
	ereturn matrix bws = `bw'
	quietly count if `touse'
	ereturn scalar N = r(N)
	ereturn scalar p = `p'
	ereturn scalar kink_unknown = `kink1'
	ereturn local cmd "rdbw2d_dist"
	ereturn local rdmodel "rdbw2d_dist"
	ereturn local kernel "`kernel'"
	ereturn local bwselect "`bwselect'"
	ereturn local bwparam "`bwparam'"
	ereturn local vce "`vce'"
	ereturn local masspoints "`masspoints'"
	ereturn local fuzzy "`fuzzy'"

	di as text _newline "Distance-based bandwidth selection"
	matlist e(bws), names(columns) format(%10.4f)
end
