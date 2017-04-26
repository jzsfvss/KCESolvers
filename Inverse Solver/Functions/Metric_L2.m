function d = Metric_L2(x, y)

%[s, n] = sumsqr(x-y);
%[s, n] = MySumSqr(x-y);
[s, n] = MySumSqr_mex(x-y);
d = sqrt(s/n);