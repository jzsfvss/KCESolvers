function [ xopt, yopt, convd ] = OptAlgGen_Master(n, f, xy, magfac, yeps, miter)

switch n
case 1
	[ xopt, yopt, convd ] = OptAlgGen_NelderMeadLog(f, xy, yeps, miter);
case 2
	% [ xopt, yopt, convd ] = OptAlgGen_CuckooSearch(f, xy, magfac, yeps, miter);
	[ xopt, yopt, convd ] = OptAlgGen_CuckooSearchLog(f, xy, magfac, yeps, miter);
case 3
	[ xopt, yopt, convd ] = OptAlgGen_FlowerPollinationLog(f, xy, magfac, yeps, miter);
case 4
	[ xopt, yopt, convd ] = OptAlgGen_PatternSearchLogML(f, xy, magfac, yeps, miter);
case 5
	[ xopt, yopt, convd ] = OptAlgGen_HarmonySearchLog(f, xy, magfac, yeps, miter);
case 6
	[ xopt, yopt, convd ] = OptAlgGen_FminconLogML(1, f, xy, magfac, yeps, miter);
case 7
	[ xopt, yopt, convd ] = OptAlgGen_FminconLogML(2, f, xy, magfac, yeps, miter);
case 8
	[ xopt, yopt, convd ] = OptAlgGen_FminconLogML(3, f, xy, magfac, yeps, miter);
case 9
	[ xopt, yopt, convd ] = OptAlgGen_ACDLog(f, xy, magfac, yeps, miter);
case 10
	[ xopt, yopt, convd ] = OptAlgGen_FminconLogML(4, f, xy, magfac, yeps, miter);
end