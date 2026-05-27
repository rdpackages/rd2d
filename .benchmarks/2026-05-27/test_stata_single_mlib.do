version 16.0
set more off
set type double

args pkg_dir
if "`pkg_dir'" == "" {
    display as error "usage: StataMP-64 /e do test_stata_single_mlib.do <package_dir>"
    exit 198
}

adopath ++ "`pkg_dir'"

capture confirm file "`pkg_dir'/rd2d_functions.do"
if !_rc {
    display as error "single-mlib test package must not include rd2d_functions.do"
    exit 198
}

import delimited using "`pkg_dir'/rd2d_data.csv", clear varnames(1)
recast double _all
capture rename x_1 x1
capture rename x_2 x2
capture rename Y y
capture rename Z_1 z1
capture rename Z_2 z2

local bpoints 0 40 0 0

macro drop RD2D_MATA_LOADED
rd2d y x1 x2 assignment, b(`bpoints') fitmethod(separate) masspoints(off)
rdbw2d y x1 x2 assignment, b(`bpoints') fitmethod(separate) masspoints(off)

generate double dist1 = sqrt((x1 - 0)^2 + (x2 - 40)^2) * (2 * assignment - 1)
generate double dist2 = sqrt((x1 - 0)^2 + (x2 - 0)^2) * (2 * assignment - 1)
rd2d_dist y dist1 dist2, b(`bpoints') fitmethod(separate) masspoints(off)
rdbw2d_dist y dist1 dist2, b(`bpoints') fitmethod(separate) masspoints(off)

mata: st_numscalar("rd2d_single_mlib_ok", rd2d_mlib_loaded())
assert rd2d_single_mlib_ok == 1
