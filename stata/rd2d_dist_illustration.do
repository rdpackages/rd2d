********************************************************************************
** RD2D Stata Package
** Distance-based empirical illustration
********************************************************************************
clear all
set more off
set linesize 100

set seed 123
set obs 800

gen x = runiform()*2 - 1
gen d = x >= 0
gen y = 1 + 2*x + 1.5*d + rnormal()*.5

********************************************************************************
** Sharp distance-based design
********************************************************************************
rdbw2d_dist y x, b(0 0) bwcheck(10)
rd2d_dist  y x, b(0 0) h(.5) bwcheck(10)

********************************************************************************
** Fuzzy distance-based design
********************************************************************************
gen takeup = runiform() < cond(d, .8, .2)
gen yf = 1 + 2*x + 1.5*takeup + rnormal()*.5

rdbw2d_dist yf x, b(0 0) fuzzy(takeup) bwcheck(10)
rd2d_dist  yf x, b(0 0) h(.5) fuzzy(takeup) bwcheck(10)

********************************************************************************
** Alias commands
********************************************************************************
rd2d_distance yf x, b(0 0) h(.5) fuzzy(takeup) bwcheck(10)
rdbw2d_distance yf x, b(0 0) fuzzy(takeup) bwcheck(10)
