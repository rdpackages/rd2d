********************************************************************************
* RD2D STATA PACKAGE -- rd2d
* Authors: Matias D. Cattaneo, Rocio Titiunik, Ruiqi Rae Yu
********************************************************************************
*!version 0.1.0  2026-05-19

capture program drop rd2d
program define rd2d, eclass
	version 16.0
	syntax varlist(min=4 max=4 numeric) [if] [in], B(numlist min=2) ///
		[ H(numlist) Fuzzy(varname numeric) DERiv(numlist max=2) ///
		  TANGvec(numlist) P(integer 1) Q(integer -1) ///
		  KERnel(string) KERNELtype(string) VCE(string) ///
		  MASSpoints(string) CLuster(varname numeric) Level(real 95) ///
		  SIDE(string) REPp(integer 1000) BWSELECT(string) BWPARAM(string) ///
		  Method(string) BWCHECK(integer -1) SCALEREGUL(real 3) ///
		  SCALEBIASCRCT(real 1) STDVARS(string) PARAMSOther(string) ///
		  PARAMSCov(string) ]

	marksample touse
	tokenize `varlist'
	local y `1'
	local x1 `2'
	local x2 `3'
	local d `4'

	markout `touse' `y' `x1' `x2' `d'
	if "`fuzzy'" != "" markout `touse' `fuzzy'
	if "`cluster'" != "" markout `touse' `cluster'

	if "`kernel'" == "" local kernel "tri"
	if "`kerneltype'" == "" local kerneltype "prod"
	if "`vce'" == "" local vce "hc1"
	if "`masspoints'" == "" local masspoints "check"
	if "`side'" == "" local side "two"
	if "`bwselect'" == "" local bwselect "mserd"
	if "`bwparam'" == "" local bwparam "main"
	if "`method'" == "" local method "dpi"
	if `bwcheck' < 0 local bwcheck = 50 + `p' + 1
	local stdflag = 1
	if lower("`stdvars'") == "off" local stdflag = 0

	if "$RD2D_MATA_LOADED" != "1" {
		local rd2d_loadonly 1
		quietly findfile rd2d_functions.do
		quietly do "`r(fn)'"
		global RD2D_MATA_LOADED 1
	}

	tempname main bw main0 main1 itt fs itt0 itt1 fs0 fs1
	tempname Vmain Vitt Vfs Vitt0 Vitt1 Vfs0 Vfs1
	mata: rd2d_location_stata("`main'", "`bw'", "`main0'", "`main1'", ///
		"`itt'", "`fs'", "`itt0'", "`itt1'", "`fs0'", "`fs1'", ///
		"`Vmain'", "`Vitt'", "`Vfs'", "`Vitt0'", "`Vitt1'", "`Vfs0'", "`Vfs1'", ///
		"`y'", "`x1'", "`x2'", "`d'", ///
		"`fuzzy'", "`cluster'", "`touse'", st_local("b"), st_local("h"), ///
		st_local("deriv"), st_local("tangvec"), `p', `q', "`kernel'", ///
		"`kerneltype'", "`vce'", `level', "`side'", "`bwselect'", ///
		`bwcheck', `stdflag')

	local maincols b1 b2 estimate_p std_err_p estimate_q std_err_q t_value p_value ci_lower ci_upper h01 h02 h11 h12 N_Co N_Tr
	local bwcols b1 b2 h01 h02 h11 h12 N_Co N_Tr
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
	ereturn local cmd "rd2d"
	if "`fuzzy'" == "" ereturn local rdmodel "rd2d"
	else ereturn local rdmodel "fuzzy rd2d"
	ereturn local kernel "`kernel'"
	ereturn local kerneltype "`kerneltype'"
	ereturn local vce "`vce'"
	ereturn local bwselect "`bwselect'"
	ereturn local bwparam "`bwparam'"
	ereturn local method "`method'"
	ereturn local masspoints "`masspoints'"
	ereturn local fuzzy "`fuzzy'"

	di as text _newline "Location-based boundary discontinuity estimates"
	matlist e(main), names(columns) format(%10.4f)
	if "`fuzzy'" != "" {
		di as text _newline "Reduced-form outcome discontinuity"
		matlist e(itt), names(columns) format(%10.4f)
		di as text _newline "First-stage treatment discontinuity"
		matlist e(fs), names(columns) format(%10.4f)
	}
end
