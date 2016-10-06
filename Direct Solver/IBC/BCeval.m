function bc = BCeval(t, mth)

switch mth
case 1
	bc = BC_MASKE(t);
case 2
	bc = BC_SimNECEEM(t);
case 3
	bc = BC_NECEEM(t);
case 4
	bc = BC_cNECEEM(t);
case 5
	bc = BC_SweepCE(t);
case 6
	bc = BC_sSweepCE(t);
case 7
	bc = BC_sSweepCEEM(t);
case 8
	bc = BC_ECEEM(t);
otherwise
	bc = BC_ppKCE(t);
end