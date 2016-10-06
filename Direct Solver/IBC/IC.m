function ic = IC(mth)

switch mth
case 1
	ic = @IC_MASKE;
case 2
	ic = @IC_SimNECEEM;
case 3
	ic = @IC_NECEEM;
case 4
	ic = @IC_cNECEEM;
case 5
	ic = @IC_SweepCE;
case 6
	ic = @IC_sSweepCE;
case 7
	ic = @IC_sSweepCEEM;
case 8
	ic = @IC_ECEEM;
case 9
	ic = @IC_ppKCE;
end