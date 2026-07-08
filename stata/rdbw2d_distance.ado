*!version 1.0.1  2026-07-08
capture program drop rdbw2d_distance
program define rdbw2d_distance, eclass
	version 16.0
	rdbw2d_dist `0'
	ereturn local cmd "rdbw2d_distance"
end
