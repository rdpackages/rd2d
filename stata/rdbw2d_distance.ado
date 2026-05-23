*!version 0.2.0  2026-05-23
capture program drop rdbw2d_distance
program define rdbw2d_distance, eclass
	version 16.0
	rdbw2d_dist `0'
	ereturn local cmd "rdbw2d_distance"
end
