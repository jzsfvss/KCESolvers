function ic = ICeval(x, n, mth)

switch mth
case 1
	ic = IC_MASKE(x, n);
case 2
	ic = IC_SimNECEEM(x, n);
case 3
	ic = IC_NECEEM(x, n);
case 4
	ic = IC_cNECEEM(x, n);
case 5
	ic = IC_SweepCE(x, n);
case 6
	ic = IC_sSweepCE(x, n);
case 7
	ic = IC_sSweepCEEM(x, n);
case 8
	ic = IC_ECEEM(x, n);
otherwise
	ic = IC_ppKCE(x, n);
end