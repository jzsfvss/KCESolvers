function [ xopt, yopt, convd ] = OptAlg_NelderMead(f, xy, yeps, miter)
% dim = 2
% y1 = f(x1) <= y2 = f(x2) <= y3 = f(x3)

%global krec

% Optimization parameters:
pref = 1; % Reflection factor.
pext = 2; % Extension factor.
pcon1 = 0.5; % Contraction factor 1.
pcon2 = 0.5; % Contraction factor 2.

% Initialization:
x = xy(1:3, 1:2);
y = xy(1:3, 3)';

% Optimization:
iter = 0;
while ((y(1) > yeps) && (iter < miter))

x0 = (x(1,:) + x(2,:))/2; % Centroid.
xr = x0 + pref*(x0 - x(3,:)); % Reflect the worst point about the centroid of the others.
yr = f(xr);
x(4,:) = xr;
y(4) = yr;
%krec = [ krec; xr, yr ];

if (yr < y(1))
	x(5,:) = x0 + pext*(xr - x0); % Extended point.
	y(5) = f(x(5,:));
	%krec = [ krec; x(5,:), y(5) ];
end

if (yr >= y(2))
	xc = x0 + pcon1*(x(3,:) - x0); % Contracting x3 about the centroid.
	yc = f(xc);
	%krec = [ krec; xc, yc ];
	if ((yr < y(3)) || (yc < y(3)))
		x(5,:) = xc;
		y(5) = yc;
	else
	for k = 2:3
		x(k,:) = x(1,:) + pcon2*(x(k,:) - x(1,:)); % Contracting x2, x3.
		y(k) = f(x(k,:));
		%krec = [ krec; x(k,:), y(k) ];
	end
	end
end

[ y, ind ] = sort(y);
x = x(ind, :);
iter = iter + 1;

end

% Output:
xopt = x(1,:);
yopt = y(1);
convd = (yopt <= yeps);