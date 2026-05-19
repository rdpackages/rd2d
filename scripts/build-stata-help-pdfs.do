version 16
clear all
set more off

args repo_root
if `"`repo_root'"' == "" {
    local repo_root "`c(pwd)'"
}
local repo_root : subinstr local repo_root "\" "/", all
local stata_dir "`repo_root'/stata"

cd "`stata_dir'"

local helpfiles rd2d rdbw2d rd2d_dist rdbw2d_dist rd2d_distance rdbw2d_distance
foreach helpfile of local helpfiles {
    display as text "Building `helpfile'.pdf"
    translate "`helpfile'.sthlp" "`helpfile'.pdf", replace translator(smcl2pdf)
}

display as text "Stata help PDF generation complete."
