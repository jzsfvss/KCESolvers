function op = Plugs()

global l;
global u0;

mth1 = 2;
disp(' ');
[ k, v, D, l, l2, u0, u022, xmax, tmax, J, I, Lini, Tini, psnm ] = ParsFile('parameters.xls', mth1, 1);

x = (-0.25*l):(l/1000):(2.25*l);
lx = length(x);
y = zeros(5,lx);

ICf = IC(1);
ord = [ 1, 5, 3, 2, 4 ];
for n = 1:5
for i = 1:lx
	ic = ICf(x(i), ord(n));
	y(n,i) = ic(1);
end
end

clf;
hold on;
%vline(0, 'k-');
%vline(l, 'k-');
col = [ 'c', 'b', 'g', 'r', 'y' ];
wid = 2*[ 2, 2, 2, 4, 2 ];
maxy = max(max(y));
axis([ x(1), x(lx), 0, 1.1*maxy ]);
for n = 1:5
	plot(x, y(n,:), '-', 'Color', col(n), 'LineWidth', wid(n));
end
%title([ 'Plug Profiles' ])
legend(DensityName(ord(1)), DensityName(ord(2)), DensityName(ord(3)), DensityName(ord(4)), DensityName(ord(5)));
xlabel('Position on the capillary (m)');
ylabel('Concentration (mol / m^3)');
hold off;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(gca, 'Color', 'None');
glg = findobj(gcf,'Type','axes','Tag','legend');
set(glg, 'Color', 'None');

disp(' ');
disp('Plotting and exporting figure to: ./Exports/plugs.png');
if (exist([ pwd, '\Exports']) == 0)
	system('mkdir Exports');
end
export_fig './Exports/plugs' -png -transparent -m3;
disp(' ');
disp('Finished.');
disp(' ');

op = 1;