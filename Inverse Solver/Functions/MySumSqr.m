function [ y, n ] = MySumSqr(x)

y = sum(sum(x.*x));
sz = size(x);
n = max(sz);