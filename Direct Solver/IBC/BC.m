function bc = BC(mth)

switch mth
case 1
	bc = @BC_MASKE;
case 2
	bc = @BC_SimNECEEM;
case 3
	bc = @BC_NECEEM;
case 4
	bc = @BC_cNECEEM;
case 5
	bc = @BC_SweepCE;
case 6
	bc = @BC_sSweepCE;
case 7
	bc = @BC_sSweepCEEM;
case 8
	bc = @BC_ECEEM;
case 9
	bc = @BC_ppKCE;
end