********************************************************************************
* RD2D STATA PACKAGE -- rd2d_dist
* Authors: Matias D. Cattaneo, Rocio Titiunik, Ruiqi Rae Yu
********************************************************************************
*!version 0.1.0  2026-05-19

capture program drop rd2d_dist
program define rd2d_dist, eclass
	version 16.0
	syntax varlist(min=2 numeric) [if] [in] ///
		[, B(numlist) H(numlist) Fuzzy(varname numeric) P(integer 1) ///
		   Q(integer -1) KINKUNKnown(string) KINKPOSition(numlist) ///
		   KERnel(string) Level(real 95) CBANDS SIDE(string) REPp(integer 1000) ///
		   BWSELECT(string) BWPARAM(string) PARAMSOther(string) PARAMSCov(string) ///
		   VCE(string) BWCHECK(integer -1) MASSpoints(string) ///
		   CLuster(varname numeric) SCALEREGUL(real 1) CQT(real .5) ]

	marksample touse
	gettoken y distvars : varlist
	markout `touse' `y' `distvars'
	if "`fuzzy'" != "" markout `touse' `fuzzy'
	if "`cluster'" != "" markout `touse' `cluster'

	if "`kernel'" == "" local kernel "tri"
	if "`side'" == "" local side "two"
	if "`bwselect'" == "" local bwselect "mserd"
	if "`bwparam'" == "" local bwparam "main"
	if "`vce'" == "" local vce "hc1"
	if "`masspoints'" == "" local masspoints "check"
	if `bwcheck' < 0 local bwcheck = 50 + `p' + 1

	local kink1 = 0
	local kink2 = 0
	local ku = lower(strtrim("`kinkunknown'"))
	if inlist("`ku'", "on", "yes", "true", "1") {
		local kink1 = 1
		local kink2 = 1
	}
	else if "`ku'" != "" & !inlist("`ku'", "off", "no", "false", "0") {
		tokenize "`ku'"
		local kink1 = real("`1'")
		local kink2 = real("`2'")
		if "`2'" == "" local kink2 = `kink1'
	}

	if "$RD2D_MATA_LOADED" != "1" {
		local rd2d_loadonly 1
		quietly findfile rd2d_functions.do
		quietly do "`r(fn)'"
		global RD2D_MATA_LOADED 1
	}

	tempname main bw main0 main1 itt fs itt0 itt1 fs0 fs1
	tempname Vmain Vitt Vfs Vitt0 Vitt1 Vfs0 Vfs1
	mata: rd2d_distance_stata("`main'", "`bw'", "`main0'", "`main1'", ///
		"`itt'", "`fs'", "`itt0'", "`itt1'", "`fs0'", "`fs1'", ///
		"`Vmain'", "`Vitt'", "`Vfs'", "`Vitt0'", "`Vitt1'", "`Vfs0'", "`Vfs1'", ///
		"`y'", "`distvars'", "`fuzzy'", ///
		"`cluster'", "`touse'", st_local("b"), st_local("h"), `p', ///
		`q', "`kernel'", "`vce'", `level', "`side'", "`bwselect'", ///
		`bwcheck', `kink1', `kink2')

	local maincols b1 b2 estimate_p std_err_p estimate_q std_err_q t_value p_value ci_lower ci_upper h0 h1 h0_rbc h1_rbc N_Co N_Tr
	local bwcols b1 b2 h0 h1 N_Co N_Tr
	capture matrix colnames `main' = `maincols'
	capture matrix colnames `bw' = `bwcols'
	capture matrix colnames `main0' = `maincols'
	capture matrix colnames `main1' = `maincols'
	capture matrix colnames `itt' = `maincols'
	capture matrix colnames `fs' = `maincols'
	capture matrix colnames `itt0' = `maincols'
	capture matrix colnames `itt1' = `maincols'
	capture matrix colnames `fs0' = `maincols'
	capture matrix colnames `fs1' = `maincols'

	ereturn clear
	ereturn matrix main = `main'
	ereturn matrix bw = `bw'
	capture confirm matrix `main0'
	if !_rc {
		if rowsof(`main0') > 0 ereturn matrix main_0 = `main0'
	}
	capture confirm matrix `main1'
	if !_rc {
		if rowsof(`main1') > 0 ereturn matrix main_1 = `main1'
	}
	capture confirm matrix `itt'
	if !_rc {
		if rowsof(`itt') > 0 ereturn matrix itt = `itt'
	}
	capture confirm matrix `fs'
	if !_rc {
		if rowsof(`fs') > 0 ereturn matrix fs = `fs'
	}
	capture confirm matrix `itt0'
	if !_rc {
		if rowsof(`itt0') > 0 ereturn matrix itt_0 = `itt0'
	}
	capture confirm matrix `itt1'
	if !_rc {
		if rowsof(`itt1') > 0 ereturn matrix itt_1 = `itt1'
	}
	capture confirm matrix `fs0'
	if !_rc {
		if rowsof(`fs0') > 0 ereturn matrix fs_0 = `fs0'
	}
	capture confirm matrix `fs1'
	if !_rc {
		if rowsof(`fs1') > 0 ereturn matrix fs_1 = `fs1'
	}
	ereturn matrix V_main = `Vmain'
	capture confirm matrix `Vitt'
	if !_rc {
		if rowsof(`Vitt') > 0 ereturn matrix V_itt = `Vitt'
	}
	capture confirm matrix `Vfs'
	if !_rc {
		if rowsof(`Vfs') > 0 ereturn matrix V_fs = `Vfs'
	}
	capture confirm matrix `Vitt0'
	if !_rc {
		if rowsof(`Vitt0') > 0 ereturn matrix V_itt_0 = `Vitt0'
	}
	capture confirm matrix `Vitt1'
	if !_rc {
		if rowsof(`Vitt1') > 0 ereturn matrix V_itt_1 = `Vitt1'
	}
	capture confirm matrix `Vfs0'
	if !_rc {
		if rowsof(`Vfs0') > 0 ereturn matrix V_fs_0 = `Vfs0'
	}
	capture confirm matrix `Vfs1'
	if !_rc {
		if rowsof(`Vfs1') > 0 ereturn matrix V_fs_1 = `Vfs1'
	}
	quietly count if `touse'
	ereturn scalar N = r(N)
	ereturn scalar p = `p'
	ereturn scalar q = cond(`q' < 0, `p' + 1, `q')
	ereturn scalar level = `level'
	ereturn scalar repp = `repp'
	ereturn scalar kink_unknown_point = `kink1'
	ereturn scalar kink_unknown_inference = `kink2'
	ereturn local cmd "rd2d_dist"
	if "`fuzzy'" == "" ereturn local rdmodel "rd2d_dist"
	else ereturn local rdmodel "fuzzy rd2d_dist"
	ereturn local kernel "`kernel'"
	ereturn local vce "`vce'"
	ereturn local bwselect "`bwselect'"
	ereturn local bwparam "`bwparam'"
	ereturn local masspoints "`masspoints'"
	ereturn local fuzzy "`fuzzy'"

	di as text _newline "Distance-based boundary discontinuity estimates"
	matlist e(main), names(columns) format(%10.4f)
	if "`fuzzy'" != "" {
		di as text _newline "Reduced-form outcome discontinuity"
		matlist e(itt), names(columns) format(%10.4f)
		di as text _newline "First-stage treatment discontinuity"
		matlist e(fs), names(columns) format(%10.4f)
	}
end
