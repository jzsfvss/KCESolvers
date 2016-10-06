function [ u, t, x, st, sx ] = DirSolverMASKE(ptp, I, J0)

% Initialization:
[ cpars, cx, cic, cbc, t, x ] = DirSolverMASKEIni_mex(ptp, I, J0);

% Solution:
u = DirSolverMASKECore_mex(cpars, cx, cic, cbc, x);
u = permute(u, [ 3, 1, 2 ]);

% Size outputs:
sz = size(u);
st = sz(2)-1;
sx = sz(3)-1;