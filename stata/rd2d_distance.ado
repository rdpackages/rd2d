*!version 0.1.0  2026-05-19
capture program drop rd2d_distance
program define rd2d_distance, eclass
	version 16.0
	rd2d_dist `0'
	ereturn local cmd "rd2d_distance"
end
