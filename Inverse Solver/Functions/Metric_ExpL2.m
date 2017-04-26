function d = Metric_ExpL2(x, y)

%[s, n] = sumsqr(x-y);
%[s, n] = MySumSqr(x-y);
[s, n] = MySumSqr_mex(x-y);
d = sqrt(s/n);
%d = exp(d) - 1;
d = (10^d) - 1;