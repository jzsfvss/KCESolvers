function nm = DensityName(n)

switch n
case 1
	nm = 'Heaviside';
case 2
	nm = 'Gaussian';
case 3
	nm = 'Shifted-Gaussian';
case 4
	nm = 'Gaussian-Heaviside';
otherwise
	nm = 'Quartic';
end