version 16.0
set more off
set type double

args repo_root mode
if "`repo_root'" == "" {
    display as error "usage: StataMP-64 /e do profile_stata.do <repo_root> [separate|joint]"
    exit 198
}
if "`mode'" == "" local mode "separate"

if "`mode'" == "separate" {
    local fitmethod "separate"
    local fm "fitmethod(separate)"
}
else if "`mode'" == "joint" {
    local fitmethod "joint"
    local fm "fitmethod(joint)"
}
else {
    display as error "unknown mode: `mode'"
    exit 198
}

adopath ++ "`repo_root'/stata"

local profile_dir "`repo_root'/.benchmarks/2026-05-27/profiles"
capture mkdir "`profile_dir'"

local timings "`profile_dir'/stata_`fitmethod'_profile_timings.csv"
global RD2D_PROFILE_TIMINGS "`timings'"
global RD2D_PROFILE_FITMETHOD "`fitmethod'"

file open tf using "`timings'", write replace
file write tf "case,fitmethod,rep,elapsed_seconds" _n
file close tf

program define append_time
    args case_name rep elapsed
    file open tf using "$RD2D_PROFILE_TIMINGS", write append
    file write tf "`case_name',$RD2D_PROFILE_FITMETHOD,`rep',`elapsed'" _n
    file close tf
end

program define time_cmd
    syntax, CASE(string) REP(integer)
    timer clear 1
    timer on 1
    quietly $RD2D_PROFILE_CMD
    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    append_time "`case'" `rep' `elapsed'
end

which rd2d
which rd2d_dist

********************************************************************************
* Large synthetic location workload
********************************************************************************
clear
set obs 4500
set seed 2026052701
generate double x1 = round(rnormal() * 1.4, .01)
generate double x2 = round(0.35 * x1 + rnormal() * 1.2, .01)
generate byte assignment = (x1 + 0.25 * x2 >= 0)
generate byte fuzzy = (runiform() < cond(assignment == 1, 0.80, 0.20))
generate int cluster = ceil(runiform() * 450)
generate double y = 1 + 0.7 * x1 - 0.4 * x2 + 0.9 * assignment + rnormal() * 0.6
generate double y_fuzzy = 1 + 0.7 * x1 - 0.4 * x2 + 1.2 * fuzzy + rnormal() * 0.6

local b_loc
forvalues j = 1/15 {
    local theta = -0.9 + (`j' - 1) * (1.8 / 14)
    local bx = 0.25 * `theta'
    local b_loc "`b_loc' `bx' `theta'"
}

local cmd `"rd2d y x1 x2 assignment, b(`b_loc') h(0.95) vce(hc0) masspoints(off) bwcheck(0) paramsother(main.0 main.1) paramscov(main main.0 main.1) `fm'"'
quietly `cmd'
global RD2D_PROFILE_CMD `"`cmd'"'
forvalues r = 1/3 {
    time_cmd, case("loc_sharp_cov") rep(`r')
}

local cmd `"rd2d y_fuzzy x1 x2 assignment, b(`b_loc') h(0.95) fuzzy(fuzzy) cluster(cluster) vce(hc1) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
quietly `cmd'
global RD2D_PROFILE_CMD `"`cmd'"'
forvalues r = 1/3 {
    time_cmd, case("loc_fuzzy_cluster_cov") rep(`r')
}

local cmd `"rdbw2d y x1 x2 assignment, b(`b_loc') vce(hc1) masspoints(off) bwcheck(0) stdvars(off) `fm'"'
quietly `cmd'
global RD2D_PROFILE_CMD `"`cmd'"'
forvalues r = 1/3 {
    time_cmd, case("loc_bw") rep(`r')
}

********************************************************************************
* Large synthetic distance workload
********************************************************************************
clear
set obs 7000
set seed 2026052702
generate double d0 = runiform() * 2.5 - 1.25
generate byte assignment = (d0 >= 0)
generate byte fuzzy = (runiform() < cond(assignment == 1, 0.80 + 0.04 * d0, 0.20 + 0.04 * d0))
generate int cluster = ceil(runiform() * 550)
generate double y = 1 + 2 * d0 + 2.4 * assignment + rnormal() * 0.3
generate double y_fuzzy = 1 + 2 * d0 + 1.5 * fuzzy + rnormal() * 0.3

local dvars
local b_dist
forvalues j = 1/15 {
    local shift = -0.3 + (`j' - 1) * (0.6 / 14)
    generate double d`j' = d0 - `shift'
    local dvars "`dvars' d`j'"
    local b_dist "`b_dist' 0 0"
}

local cmd `"rd2d_dist y `dvars', b(`b_dist') h(0.5) p(1) q(2) vce(hc0) kernel(tri) masspoints(off) bwcheck(0) paramsother(main.0 main.1) paramscov(main main.0 main.1) `fm'"'
quietly `cmd'
global RD2D_PROFILE_CMD `"`cmd'"'
forvalues r = 1/3 {
    time_cmd, case("dist_sharp_cov") rep(`r')
}

local cmd `"rd2d_dist y_fuzzy `dvars', b(`b_dist') h(0.5) p(1) q(2) fuzzy(fuzzy) cluster(cluster) vce(hc1) kernel(tri) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
quietly `cmd'
global RD2D_PROFILE_CMD `"`cmd'"'
forvalues r = 1/3 {
    time_cmd, case("dist_fuzzy_cluster_cov") rep(`r')
}

local cmd `"rdbw2d_dist y `dvars', b(`b_dist') p(1) vce(hc1) kernel(tri) masspoints(off) bwcheck(0) `fm'"'
quietly `cmd'
global RD2D_PROFILE_CMD `"`cmd'"'
forvalues r = 1/3 {
    time_cmd, case("dist_bw") rep(`r')
}

exit
