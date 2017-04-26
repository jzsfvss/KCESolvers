function d = Metric_L2W(x, y, w)

%[s, n] = sumsqr(x-y);
%[s, n] = MySumSqr(x-y);
[s, n] = MySumSqr_mex(x-y);
d = sqrt(sum(w.*s)/sum(w));