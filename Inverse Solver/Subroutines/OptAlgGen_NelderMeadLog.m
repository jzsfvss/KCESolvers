function [ xopt, yopt, convd ] = OptAlgGen_NelderMeadLog(f0, xy, yeps, miter)
% y1 = f(x1) <= y2 = f(x2) <= ...

% Optimization parameters:
pref = 1; % Reflection factor.
pext = 2; % Extension factor.
pcon1 = 0.5; % Contraction factor 1.
pcon2 = 0.5; % Contraction factor 2.

% Initialization:
f = @(x) f0(10.^x);
sz = size(xy);
d = sz(2)-1;
x = xy(1:(d+1), 1:d);
y = xy(1:(d+1), d+1)';
%zrs = 0*(1:d);
ons = 0*(1:d);
x = [ log10(x); -Inf*ons; -Inf*ons ];
y = [ y, Inf, Inf ];

% Optimization:
iter = 0;
while ((y(1) > yeps) && (iter < miter))

x0 = sum(x)/d; % Centroid.
xr = x0 + pref*(x0 - x(d+1, :)); % Reflect the worst point about the centroid of the others.
yr = f(xr);
x(d+2, :) = xr;
y(d+2) = yr;
%krec = [ krec; xr, yr ];

if (yr < y(1))
	x(d+3, :) = x0 + pext*(xr - x0); % Extended point.
	y(d+3) = f(x(d+3, :));
end

if (yr >= y(2))
	xc = x0 + pcon1*(x(d+1, :) - x0); % Contracting x3 about the centroid.
	yc = f(xc);
	if ((yr < y(d+1)) || (yc < y(d+1)))
		x(d+3, :) = xc;
		y(d+3) = yc;
	else
	for k = 2:(d+1)
		x(k, :) = x(1, :) + pcon2*(x(k, :) - x(1, :)); % Contracting x2, x3.
		y(k) = f(x(k, :));
	end
	end
end

[ y, ind ] = sort(y);
x = x(ind, :);
iter = iter + 1;

end

% Output:
xopt = 10.^x(1,:);
yopt = y(1);
convd = (yopt <= yeps);