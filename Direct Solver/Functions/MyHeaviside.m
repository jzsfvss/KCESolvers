function y = MyHeaviside(x)

global xmax

% Factors:
meps = 0.05; % Width of interpolation from y = 0 (truncation).

% Initialization:
y = 0*x;

% Markers:
m0 = 0;
m1 = m0 + meps; % Spline interpolation on the [ m0, m1 ] interval.

% Indices:
I1 = find((x >= m0) & (x < m1)); % Interpolated left-truncation.
I2 = find(x >= m1);

% Evaluation of the Heaviside part:
if (~isempty(I2))
	y(I2) = 1;
end

% Evaluation of the spline part:
if (~isempty(I1))

x1 = m0;
x2 = m1;
y1 = 0;
y2 = 1;
sl1 = 0;
sl2 = 0;
% y(I1) = MyHermite(x(I1), x1, y1, sl1, x2, y2, sl2);
y(I1) = y1 + ((y2 - y1)/(x2 - x1))*(x(I1) - x1); % Linear interpolation.

end