function [ u, t, x, st, sx ] = DirSolver(mth1, mth2, ptp, I, J)

switch (mth1)

case 1 % Standard KCE:

[ u, t, x, st, sx ] = DirSolverKCE(mth2 + 2, ptp, I, J);

case 2 % MASKE:

switch (mth2)
case 1 % Exact:
	[ u, t, x, st, sx ] = ExactSolMASKE(IC(1), ptp, I, 1);
case 2 % Numerical:
	[ u, t, x, st, sx ] = DirSolverMASKE(ptp, I, J);
end

case 3 % Simplified NECEEM:

switch (mth2)
case 1 % Exact:
	[ u, t, x, st, sx ] = ExactSolSimNECEEM(IC(2), ptp, I, 1, 0.01);
case 2 % Numerical:
	[ u, t, x, st, sx ] = DirSolverKCE(2, ptp, I, J);
end

end