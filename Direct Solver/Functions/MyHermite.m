function y = MyHermite(x, x1, y1, m1, x2, y2, m2)
% Cubic Hermite spline interpolation.

% Rescaling:
t = (x - x1)/(x2 - x1);

% Hermite basis functions:
h00t = (1 + 2*t).*(1 - t).*(1 - t);
h10t = t.*(1 - t).*(1 - t);
h01t = t.*t.*(3 - 2*t);
h11t = t.*t.*(t - 1);

% Interpolated values:
y = h00t*y1 + h01t*y2 + (h10t*m1 + h11t*m2)*(x2-x1);