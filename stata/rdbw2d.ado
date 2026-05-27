********************************************************************************
* RD2D STATA PACKAGE -- rdbw2d
* Authors: Matias D. Cattaneo, Rocio Titiunik, Ruiqi Rae Yu
********************************************************************************
*!version 1.0.0  2026-05-26

capture program drop rdbw2d
program define rdbw2d, eclass
	version 16.0
	syntax varlist(min=4 max=4 numeric) [if] [in], B(string asis) ///
		[ Fuzzy(varname numeric) DERiv(numlist max=2) TANGvec(numlist) ///
		  P(integer 1) KERnel(string) KERNELtype(string) ///
		  BWSELECT(string) BWPARAM(string) Method(string) VCE(string) ///
		  BWCHECK(integer 20) MASSpoints(string) CLuster(varname numeric) ///
		  SCALEREGUL(real 1) SCALEBIASCRCT(real 1) STDVARS(string) ///
		  FITMethod(string) COVSEFF(varlist numeric) ]

	marksample touse
	tokenize `varlist'
	local y `1'
	local x1 `2'
	local x2 `3'
	local d `4'

	markout `touse' `y' `x1' `x2' `d'
	if "`fuzzy'" != "" markout `touse' `fuzzy'
	if "`covseff'" != "" markout `touse' `covseff'
	if "`cluster'" != "" markout `touse' `cluster'

	if "`kernel'" == "" local kernel "tri"
	if "`kerneltype'" == "" local kerneltype "prod"
	if "`bwselect'" == "" local bwselect "mserd"
	if "`bwparam'" == "" local bwparam "main"
	if "`method'" == "" local method "dpi"
	if "`vce'" == "" local vce "hc1"
	if "`masspoints'" == "" local masspoints "check"
	if "`fitmethod'" == "" local fitmethod "joint"
	local fitmethod = lower("`fitmethod'")
	if !inlist("`fitmethod'", "joint", "separate") {
		di as error "fitmethod() must be joint or separate"
		exit 198
	}
	local stdflag = 1
	if lower("`stdvars'") == "off" local stdflag = 0

	if "$RD2D_MATA_LOADED" != "1" {
		tempname rd2d_mlib_ok
		capture quietly mata: mata mlib index
		capture quietly mata: st_numscalar("`rd2d_mlib_ok'", rd2d_mlib_loaded())
		if _rc {
			local rd2d_mlib_rc = _rc
			di as error "rd2d Mata library lrd2d.mlib was not found or could not be loaded"
			di as error "reinstall rd2d or rebuild lrd2d.mlib from rd2d_functions.do"
			exit `rd2d_mlib_rc'
		}
		global RD2D_MATA_LOADED 1
	}

	tempname bw
	mata: rdbw2d_location_stata("`bw'", "`y'", "`x1'", "`x2'", "`d'", ///
		"`fuzzy'", "`cluster'", "`touse'", st_local("b"), ///
		st_local("deriv"), st_local("tangvec"), `p', ///
		"`kernel'", "`kerneltype'", "`bwselect'", "`bwparam'", ///
		"`method'", "`vce'", `bwcheck', "`masspoints'", ///
		`scaleregul', `scalebiascrct', `stdflag', "`fitmethod'", "`covseff'")

	matrix colnames `bw' = b1 b2 h01 h02 h11 h12 N_Co N_Tr
	ereturn clear
	ereturn matrix bws = `bw'
	quietly count if `touse'
	ereturn scalar N = r(N)
	ereturn scalar p = `p'
	ereturn local cmd "rdbw2d"
	ereturn local rdmodel "rdbw2d"
	ereturn local kernel "`kernel'"
	ereturn local kerneltype "`kerneltype'"
	ereturn local bwselect "`bwselect'"
	ereturn local bwparam "`bwparam'"
	ereturn local method "`method'"
	ereturn local vce "`vce'"
	ereturn local fitmethod "`fitmethod'"
	ereturn local masspoints "`masspoints'"
	ereturn local fuzzy "`fuzzy'"
	ereturn local covseff "`covseff'"

	di as text _newline "Location-based bandwidth selection"
	matlist e(bws), names(columns) format(%10.4f)
end
