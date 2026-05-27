version 16.0
set more off
set type double

args repo_root mode
if "`repo_root'" == "" {
    display as error "usage: StataMP-64 /e do bench_stata.do <repo_root> <mode>"
    exit 198
}
if "`mode'" == "" local mode "public"

local bench_dir "`repo_root'/.benchmarks/2026-05-27"
local input_dir "`bench_dir'/inputs"
local result_dir "`bench_dir'/results"
capture mkdir "`result_dir'"

if "`mode'" == "public" {
    local source "public"
    local version_label "0.2.0"
    local fm
    local fitmethod "implicit_separate"
}
else if "`mode'" == "local_separate" {
    adopath ++ "`repo_root'/stata"
    local source "local_separate"
    local version_label "1.0.0"
    local fm "fitmethod(separate)"
    local fitmethod "separate"
}
else if "`mode'" == "local_joint" {
    adopath ++ "`repo_root'/stata"
    local source "local_joint"
    local version_label "1.0.0"
    local fm "fitmethod(joint)"
    local fitmethod "joint"
}
else {
    display as error "unknown mode: `mode'"
    exit 198
}

which rd2d
which rd2d_dist

local results "`result_dir'/stata_`mode'_results.csv"
local timings "`result_dir'/stata_`mode'_timings.csv"
local metadata "`result_dir'/stata_`mode'_metadata.csv"
global RD2D_BENCH_RESULTS "`results'"
global RD2D_BENCH_TIMINGS "`timings'"
global RD2D_BENCH_SOURCE "`source'"
global RD2D_BENCH_VERSION "`version_label'"
global RD2D_BENCH_FITMETHOD "`fitmethod'"

file open rf using "`results'", write replace
file write rf "platform,source,version,fitmethod,case,call,output,row,column,value" _n
file close rf

file open tf using "`timings'", write replace
file write tf "platform,source,version,fitmethod,case,rep,elapsed_seconds" _n
file close tf

file open mf using "`metadata'", write replace
file write mf "platform,source,mode,version,fitmethod" _n
file write mf "Stata,`source',`mode',`version_label',`fitmethod'" _n
file close mf

mata:
void rd2d_bench_append_matrix(
    string scalar file,
    string scalar source,
    string scalar version_label,
    string scalar fitmethod,
    string scalar case_name,
    string scalar call_name,
    string scalar output_name,
    string scalar matrix_name)
{
    real matrix M
    string matrix cs, rs
    real scalar fh, i, j
    M = st_matrix(matrix_name)
    if (rows(M) == 0 | cols(M) == 0) return
    cs = st_matrixcolstripe(matrix_name)
    rs = st_matrixrowstripe(matrix_name)
    fh = fopen(file, "a")
    for (i = 1; i <= rows(M); i++) {
        for (j = 1; j <= cols(M); j++) {
            fput(fh, "Stata," + source + "," + version_label + "," + fitmethod + "," +
                case_name + "," + call_name + "," + output_name + "," +
                (cols(rs) >= 2 ? rs[i,2] : strofreal(i)) + "," +
                (cols(cs) >= 2 ? cs[j,2] : strofreal(j)) + "," +
                sprintf("%21.15g", M[i,j]))
        }
    }
    fclose(fh)
}
end

program define append_mats
    args case_name
    local call_name "main_call"
    foreach item in main bw bws itt fs main_0 main_1 itt_0 itt_1 fs_0 fs_1 V_main V_itt V_fs V_main_0 V_main_1 V_itt_0 V_itt_1 V_fs_0 V_fs_1 {
        capture confirm matrix e(`item')
        if !_rc {
            local outname "`item'"
            local outname : subinstr local outname "_" ".", all
            if substr("`item'", 1, 2) == "V_" {
                local outname "cov.`=substr("`outname'", 3, .)'"
            }
            mata: rd2d_bench_append_matrix("$RD2D_BENCH_RESULTS", "$RD2D_BENCH_SOURCE", "$RD2D_BENCH_VERSION", "$RD2D_BENCH_FITMETHOD", "`case_name'", "`call_name'", "`outname'", "e(`item')")
        }
    }
end

program define append_time
    args case_name rep elapsed
    file open tf using "$RD2D_BENCH_TIMINGS", write append
    file write tf "Stata,$RD2D_BENCH_SOURCE,$RD2D_BENCH_VERSION,$RD2D_BENCH_FITMETHOD,`case_name',`rep',`elapsed'" _n
    file close tf
end

program define time_cmd
    syntax, CASE(string) REP(integer)
    timer clear 1
    timer on 1
    quietly $RD2D_BENCH_CMD
    timer off 1
    quietly timer list 1
    local elapsed = r(t1)
    append_time "`case'" `rep' `elapsed'
end

local b_synth_loc "-0.2 -0.3 0 0 0.3 0.2"
local b_synth_dist "0 0"

import delimited using "`input_dir'/jasa_location_eval.csv", clear varnames(1) encoding("utf-8")
local b_jasa
forvalues i = 1/`=_N' {
    local b_jasa `b_jasa' `=x1[`i']' `=x2[`i']'
}

import delimited using "`input_dir'/jss_location_eval.csv", clear varnames(1) encoding("utf-8")
local b_jss_loc
forvalues i = 1/`=_N' {
    local b_jss_loc `b_jss_loc' `=x1[`i']' `=x2[`i']'
}

import delimited using "`input_dir'/joe_distance_eval.csv", clear varnames(1) encoding("utf-8")
local b_joe_dist
forvalues i = 1/`=_N' {
    local b_joe_dist `b_joe_dist' `=x1[`i']' `=x2[`i']'
}

import delimited using "`input_dir'/jss_distance_eval.csv", clear varnames(1) encoding("utf-8")
local b_jss_dist
forvalues i = 1/`=_N' {
    local b_jss_dist `b_jss_dist' `=x1[`i']' `=x2[`i']'
}

********************************************************************************
* Synthetic location
********************************************************************************
import delimited using "`input_dir'/synthetic_location.csv", clear varnames(1)
local cmd `"rd2d y x1 x2 assignment, b(`b_synth_loc') h(0.95) vce(hc0) masspoints(off) bwcheck(0) paramsother(main.0 main.1) paramscov(main main.0 main.1) `fm'"'
quietly `cmd'
append_mats "synth_loc_sharp_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("synth_loc_sharp_fixed") rep(`r')
}

local cmd `"rd2d y_fuzzy x1 x2 assignment, b(`b_synth_loc') h(0.95) fuzzy(fuzzy) vce(hc0) masspoints(off) bwcheck(0) paramsother(itt.0 itt.1 fs.0 fs.1) paramscov(main itt fs itt.0 itt.1 fs.0 fs.1) `fm'"'
quietly `cmd'
append_mats "synth_loc_fuzzy_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("synth_loc_fuzzy_fixed") rep(`r')
}

local cmd `"rd2d y x1 x2 assignment, b(`b_synth_loc') h(0.95) cluster(cluster) vce(hc1) masspoints(off) bwcheck(0) paramscov(main) `fm'"'
quietly `cmd'
append_mats "synth_loc_cluster_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("synth_loc_cluster_fixed") rep(`r')
}

local cmd `"rdbw2d y x1 x2 assignment, b(`b_synth_loc') vce(hc1) masspoints(off) bwcheck(0) stdvars(off) `fm'"'
quietly `cmd'
append_mats "synth_loc_bw"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("synth_loc_bw") rep(`r')
}

if "`mode'" != "public" {
    local cmd `"rd2d y_fuzzy x1 x2 assignment, b(`b_synth_loc') h(0.95) fuzzy(fuzzy) covseff(z1 z2) vce(hc0) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
    quietly `cmd'
    append_mats "synth_loc_fuzzy_covs_fixed"
    global RD2D_BENCH_CMD `"`cmd'"'
    forvalues r = 1/4 {
        time_cmd, case("synth_loc_fuzzy_covs_fixed") rep(`r')
    }
}

********************************************************************************
* Synthetic distance
********************************************************************************
import delimited using "`input_dir'/synthetic_distance.csv", clear varnames(1)
local cmd `"rd2d_dist y d1, b(`b_synth_dist') h(0.5) p(1) q(2) cbands vce(hc0) kernel(tri) masspoints(off) bwcheck(0) paramsother(main.0 main.1) paramscov(main main.0 main.1) `fm'"'
quietly `cmd'
append_mats "synth_dist_sharp_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("synth_dist_sharp_fixed") rep(`r')
}

local cmd `"rd2d_dist y_fuzzy d1, b(`b_synth_dist') h(0.5) p(1) q(2) fuzzy(fuzzy) vce(hc0) kernel(tri) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
quietly `cmd'
append_mats "synth_dist_fuzzy_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("synth_dist_fuzzy_fixed") rep(`r')
}

local cmd `"rd2d_dist y d1, b(`b_synth_dist') h(0.5) p(1) q(2) cbands cluster(cluster) vce(hc1) kernel(tri) masspoints(off) bwcheck(0) paramscov(main) `fm'"'
quietly `cmd'
append_mats "synth_dist_cluster_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("synth_dist_cluster_fixed") rep(`r')
}

local cmd `"rdbw2d_dist y d1, b(`b_synth_dist') p(1) vce(hc1) kernel(tri) masspoints(off) bwcheck(0) `fm'"'
quietly `cmd'
append_mats "synth_dist_bw"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("synth_dist_bw") rep(`r')
}

if "`mode'" != "public" {
    local cmd `"rd2d_dist y_fuzzy d1, b(`b_synth_dist') h(0.5) p(1) q(2) fuzzy(fuzzy) covseff(z1 z2) vce(hc0) kernel(tri) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
    quietly `cmd'
    append_mats "synth_dist_fuzzy_covs_fixed"
    global RD2D_BENCH_CMD `"`cmd'"'
    forvalues r = 1/4 {
        time_cmd, case("synth_dist_fuzzy_covs_fixed") rep(`r')
    }
}

********************************************************************************
* Empirical anchors
********************************************************************************
import delimited using "`input_dir'/jasa_location.csv", clear varnames(1)
local cmd `"rd2d y x1 x2 assignment, b(`b_jasa') h(20) fuzzy(fuzzy) vce(hc0) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
quietly `cmd'
append_mats "jasa_loc_fuzzy_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("jasa_loc_fuzzy_fixed") rep(`r')
}

import delimited using "`input_dir'/jss_location.csv", clear varnames(1)
local cmd `"rd2d y x1 x2 assignment, b(`b_jss_loc') h(20) fuzzy(fuzzy) vce(hc0) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
quietly `cmd'
append_mats "jss_loc_fuzzy_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("jss_loc_fuzzy_fixed") rep(`r')
}

import delimited using "`input_dir'/joe_distance.csv", clear varnames(1)
local cmd `"rd2d_dist y d1 d2 d3, b(`b_joe_dist') h(20) p(1) q(2) fuzzy(fuzzy) vce(hc0) kernel(tri) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
quietly `cmd'
append_mats "joe_dist_fuzzy_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("joe_dist_fuzzy_fixed") rep(`r')
}

import delimited using "`input_dir'/jss_distance.csv", clear varnames(1)
local cmd `"rd2d_dist y d1 d2 d3, b(`b_jss_dist') h(20) p(1) q(2) fuzzy(fuzzy) vce(hc0) kernel(tri) masspoints(off) bwcheck(0) paramscov(main itt fs) `fm'"'
quietly `cmd'
append_mats "jss_dist_fuzzy_fixed"
global RD2D_BENCH_CMD `"`cmd'"'
forvalues r = 1/4 {
    time_cmd, case("jss_dist_fuzzy_fixed") rep(`r')
}

exit
