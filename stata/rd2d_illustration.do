********************************************************************************
** RD2D Stata Package
** Location-based empirical illustration
********************************************************************************
clear all
set more off
set linesize 100

set seed 123
set obs 800

gen x1 = rnormal()
gen x2 = rnormal()
gen d  = x1 >= 0
gen y  = 3 + 2*x1 + 1.5*x2 + d + rnormal()

********************************************************************************
** Sharp location-based design
********************************************************************************
rdbw2d y x1 x2 d, b(0 0 0 1) bwcheck(10)
rd2d  y x1 x2 d, b(0 0 0 1) h(.9) masspoints(off) bwcheck(10)

********************************************************************************
** Fuzzy location-based design
********************************************************************************
gen takeup = runiform() < cond(d, .8, .2)
gen yf = 3 + 2*x1 + 1.5*x2 + 1.5*takeup + rnormal()

rdbw2d yf x1 x2 d, b(0 0 0 1) fuzzy(takeup) bwcheck(10)
rd2d  yf x1 x2 d, b(0 0 0 1) h(.9) fuzzy(takeup) bwcheck(10)

********************************************************************************
** Directional derivative example
********************************************************************************
rd2d y x1 x2 d, b(0 0) h(.9) tangvec(1 0) p(1) q(1) bwcheck(10)
