********************************************************************************
** RD2D Stata Package
** Numerical Illustration
********************************************************************************

clear all
set more off
set linesize 120
set type double
version 16.0

adopath ++ "."

* Load simulated data.
import delimited using "rd2d_data.csv", clear varnames(1)
recast double _all
capture rename x_1 x1
capture rename x_2 x2
capture rename Y y

* Set illustration inputs.
local neval = 40
if "$RD2D_ILLUSTRATION_NEVAL" != "" local neval = real("$RD2D_ILLUSTRATION_NEVAL")

* Generate boundary evaluation points.
local nhalf = ceil(`neval' / 2)
local bpoints
forvalues j = 1/`nhalf' {
	local b1 = 0
	local b2 = 40 - (`j' - 1) * 40 / `nhalf'
	local bpoints `bpoints' `b1' `b2'
}
local second = `neval' - `nhalf'
forvalues j = 1/`second' {
	local b1 = (`j' - 1) * 56 / `nhalf'
	local b2 = 0
	local bpoints `bpoints' `b1' `b2'
}

* Generate signed distances to each boundary evaluation point.
forvalues j = 1/`neval' {
	if `j' <= `nhalf' {
		local b1 = 0
		local b2 = 40 - (`j' - 1) * 40 / `nhalf'
	}
	else {
		local k = `j' - `nhalf'
		local b1 = (`k' - 1) * 56 / `nhalf'
		local b2 = 0
	}
	generate double dist`j' = sqrt((x1 - `b1')^2 + (x2 - `b2')^2) * (2 * assignment - 1)
}

* Location-based bandwidth selection.
rdbw2d y x1 x2 assignment, b(`bpoints') masspoints(off)

* Location-based fuzzy estimation.
rd2d y x1 x2 assignment, b(`bpoints') fuzzy(fuzzy) masspoints(off)

* Distance-based bandwidth selection.
rdbw2d_dist y dist1-dist`neval', b(`bpoints') masspoints(off)

* Distance-based sharp estimation.
rd2d_dist y dist1-dist`neval', b(`bpoints') masspoints(off)

* Distance-based fuzzy bandwidth selection.
rdbw2d_dist y dist1-dist`neval', b(`bpoints') fuzzy(fuzzy) bwparam(itt) masspoints(off)

* Distance-based fuzzy estimation.
rd2d_dist y dist1-dist`neval', b(`bpoints') fuzzy(fuzzy) bwparam(itt) masspoints(off)
