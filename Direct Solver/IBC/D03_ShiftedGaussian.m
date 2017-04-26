function den = D03_ShiftedGaussian(x)

%den = normpdf(x, 1.5, 0.4);
den = GaussPlugDS(x, 1.5, [ 0.4, 0.4 ]);