function den = Density(x, n)

switch n
case 1
	den = D01_Heaviside(x);
case 2
	den = D02_Gaussian(x);
case 3
	den = D03_ShiftedGaussian(x);
case 4
	den = D04_GaussianHeaviside(x);
case 5
	den = D05_Quartic(x);
otherwise
	den = D07_TrGaussian(x);
end