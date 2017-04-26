function y = GaussPlugDS(x, c, s)
% Computes a Gaussian left-truncated at c - 3.5*s or 0.

global xmax

% Factors:
% stfac = 3.5; % St. devs. factor.
meps = 0.1*s(1); % Width of interpolation from y = 0 (truncation).

% Initialization:
y = 0*x;
s1 = s(1);
s2 = s(2);

% Markers:
% m0 = max(0, c - stfac*s1);
m0 = 0;
m1 = m0 + meps; % Spline interpolation on the [ m0, m1 ] interval.

% Indices:
I1 = find((x >= m0) & (x < m1)); % Interpolated left-truncation.
I2 = find((x >= m1) & (x < c)); % Gaussian with st. dev. s1.
I3 = find(x >= c); % Gaussian with st. dev. s2.

% Evaluation of the Gaussian parts:
h1 = 1/(sqrt(2*pi)*s1);
h2 = 1/(sqrt(2*pi)*s2);
ha = 1/(sqrt(2*pi)*(s1 + s2)/2);
if (~isempty(I2))
	y(I2) = (ha/h1)*normpdf(x(I2), c, s1);
end
if (~isempty(I3))
	y(I3) = (ha/h2)*normpdf(x(I3), c, s2);
end

% Evaluation of the spline part:
if (~isempty(I1))

x1 = m0;
x2 = m1;
y1 = 0;
y2 = (ha/h1)*normpdf(x2, c, s1);
sl1 = 0;
sl2 = -x2*(ha/h1)*normpdf(x2, c, s1); % Der. at x2.
y(I1) = MyHermite(x(I1), x1, y1, sl1, x2, y2, sl2);

end