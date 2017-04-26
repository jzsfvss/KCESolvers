function LpC = SolLpC(mth1, mth2, ptp, I, J)

[ u, t, x, st, sx ] = DirSolver(mth1, mth2, ptp, I, J);
sz = size(u);
indet = sz(3);
LpC = zeros(1, sz(2));

if (mth1 ~= 2)
	Cind = 3;
else
	Cind = 2;
end

LpC(:) = max(0, u(1, :, indet) + u(Cind, :, indet));