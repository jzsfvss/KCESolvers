function plt = ErrorSurfPlot(ncntrs, stax, pdg, kon0, koff0, kEnum)

global krec
global optstinds

% Initialization:
mrsz = 10; % Plot marker size.
%pdg = 0.05; % Padding %.
%pdg = 1.5;
sz = size(krec);
n = sz(1);
%{
x = krec(:, 1);
y = krec(:, 2);
z = krec(:, 3);
%}
x = real(log10(krec(:, 1)));
y = real(log10(krec(:, 2)));
x0 = real(log10(kon0));
y0 = real(log10(koff0));
z = krec(:, 3);
%z = real(log10(krec(:, 3)));

% Triangulation:
deltri = delaunay(x, y);

% Plotting:
clf;
hold on;

[ c, h, plt ] = tricontour(deltri, x, y, -z, ncntrs);
%clabel(c, h);

if (plt == 0)

clf;
hold off;

else

% Test points:
plot(log10(krec(1:kEnum, 1)), log10(krec(1:kEnum, 2)), 'k.', 'Markersize', 1);

% Initial estimate:
plot(x0, y0, 'ko');
plot(x0, y0, 'k.', 'Markersize', mrsz);

% Optimization iterates:
loi = length(optstinds);
optinds = [ optstinds, n+1 ];
for i = 1:loi
	optindsi = optinds(i):(optinds(i+1) - 1);
	plot([ x0, x(1), x(optindsi)' ], [ y0, y(1), y(optindsi)' ], 'k-');
	plot([ x(1), x(optindsi)' ], [ y(1), y(optindsi)' ], 'k.', 'Markersize', mrsz);
end

% Extreme point:
[ zmin, imin ] = min(z);
plot(x(imin), y(imin), 'ro');
plot(x(imin), y(imin), 'r.', 'Markersize', mrsz);
lh = hline(y(imin), 'r:');
lv = vline(x(imin), 'r:');

% Axes:
if (stax == 1)
	xopt = x(optstinds(1):n);
	yopt = y(optstinds(1):n);
	dx = max(xopt) - min(xopt);
	x1 = min(xopt) - pdg*dx;
	x2 = max(xopt) + pdg*dx;
	dy = max(yopt) - min(yopt);
	y1 = min(yopt) - pdg*dy;
	y2 = max(yopt) + pdg*dy;
	axis([ x1, x2, y1, y2 ]);
end

% Labels:
xlabel('log_{10}(k_{on})');
ylabel('log_{10}(k_{off})');
th = title([ 'Contour plot of the error surface from ', num2str(n), ' evaluations with log_{10}(k_{opt}) = (', num2str(x(imin), '%.2f'), ', ', num2str(y(imin), '%.2f'), ')' ]);
%P = get(th, 'Position');
%set(th, 'Position', P.*[ 1, 0.985, 1 ]);

% K_d diagonal:
lgKd = y(imin) - x(imin);
xls = xlim;
yls = xls + lgKd;
plot(xls, yls, 'r:');

hold off;

set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(gca, 'Color', 'None');
set(gca, 'TickDir', 'out');

end