function res = InitializerTest()
%__________________________________________________________________________________________
% Global parameters
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
global k % k = (k_on, k_off) vector (m^3/mol, none).
global kppars % k and additional plug optimization parameters.
global kpparsest % Initial estimate of kppars.
global krec % Matrix recording rows of (k_on, k_off, E(k_on, k_off)).
global kpparsrec % Matrix recording kppars.

global v % v = (v_L, v_T, v_C) mobilities (m^2/V*sec).
global D % D = (D_L, D_T, D_C) diffusion coefficients (m^2/sec).
global l % l = length of the injected plug (m).
global l2 % l2 = length of the second plug (m).
global l0 % Original l.
global u0 % u0 = (L0eq, T0eq, C0eq) initial equilibrium concentrations (mol/m^3).
global u022 % u022 = T0eq_2 initial concentration (mol/m^3).
global xmax % xmax = total capillary length (m).
global tmax % tmax = total time of run (sec).

global Lini % Initial pre-equilibrium concentration (M).
global Tini % Initial pre-equilibrium concentration (M).
global mLpC0 % Magnitude of the experimental signal.
global optstinds % Starting index of each opt. alg. run in krec.
global ioptvec % Vector of inverter algorithm numbers.
global ioptvecmlt % Vector of inverter algorithm multiplicities.
global iopttime % Default runtime for each inverter algorithm.
global rtrat % Runtime ratio relative to the computer where this package was implemented.
global rdirsol % Direct solver runtime.
global numtdiv % Number of direct solver runs per optimization iteration.

global rI % All relevant indices.
global rIlp % Left peak indices (C).
global rIbr % Bridge between the peaks.
global rIrp % Right peak indices (L).
%__________________________________________________________________________________________
% Initialization
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%disp(' ');
%disp('Parameter Estimation Test for Kinetic Capillary Electrophoresis');

ptp = 6;
mth1 = 1;
mth2 = 1;
mopt = 1;
Metric = MetricFun(mopt);
iopt = 2;
krec = [];
mtht = 1;

disp(' ');
estopt = OptSel({'Estimate', 'k', 'AG plug parameters'});
disp(' ');
disp('Select a center parameter set:');
disp(' ');
[ k, v, D, l, l2, u0, u022, xmax, tmax, J, I, Lini, Tini, psnm ] = ParsFile('../Direct Solver/parameters.xls', mth1, 2);
kon0 = k(1);
koff0 = k(2);
l00 = l;

N = 1:3;
kppars0 = [ k(1), k(2) ]; % Shifted-Gaussian parameters.
u00 = u0;
for n = N % Initializing kppars for L, T, C.
	kppars0 = [ kppars0, l*1.5, l*0.4, l*0.4, l*u0(n)/(sqrt(2*pi)*(l*0.4)) ];
end
kppars = kppars0;
l0 = l;
l = 1;
l2 = 1;
u0(N) = 0*u0(N) + 1;

if (I == 0)
	disp(' ');
	I = input('Mesh divisions in t = ');
	J = 1;
end

LpC0 = SolLpC(mth1, mth2, ptp, I, J);
[ rI, rIlp, rIbr, rIrp ] = RelInd(LpC0, 0.01, ptp, iopt);
%__________________________________________________________________________________________
% Initial estimate
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
l1 = l;
l = l0;
k0 = k;
Kd = k(2)/k(1);
kmtx = [];
re = zeros(1, 4);
lkp = length(kppars);

tstnum = input('Number of test points = ');
magfac = input('Order of magnitude spread = ');

t1 = tic;
for i = 1:tstnum

pow = magfac*2*(rand(1) - 0.5);
k(1) = k0(1)*(10^pow);
k(2) = Kd*k(1);
kon0 = k(1);
koff0 = k(2);

LpC0 = SolLpC(mth1, mth2, ptp, I, J);


if (estopt == 1)

[ kon1, koff1, kppars1 ] = Initializer(LpC0, rI, ptp, mth1, mth2, []);
re(1) = 100*abs(kon1 - kon0)/kon0;
re(2) = 100*abs(koff1 - koff0)/koff0;

[ kon1a, koff1a ] = InitializerArea(LpC0);
re(3) = 100*abs(kon1a - kon0)/kon0;
re(4) = 100*abs(koff1a - koff0)/koff0;

else

[ kon1, koff1, kppars1 ] = Initializer(LpC0, rI, ptp, mth1, mth2, k);
re = 100*abs(kppars1 - kppars0)./kppars0;
re = re(3:lkp);

end

kmtx = [ kmtx; k(1), re ];

if (i == 5)
	t2 = toc(t1)/5;
	rtm = t2*(tstnum - 5);
	disp(' ');
	disp([ 'Estimated runtime = ', num2str(ceil(rtm)), ' secs = ', num2str(MyRound(rtm/60, 2)), ' mins.' ]);
end % if

end % for

l = l1;
%__________________________________________________________________________________________
% Plotting
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
[ kmtx1, ind ] = sort(kmtx(:, 1));

if (estopt == 1)

clf;
hold on;

colr = [1 0 0]; % red
coll = [0 0 1]; % blue
colg = [0 1 0]; % green
colb = [0 0 0]; % black
lnw1 = 3;
lnw2 = 1.5;
mrsz = 10;

axis([ log10(k0(1)*(10^(-magfac))), log10(k0(1)*(10^magfac)), 0, 100 ]);

plot(log10(kmtx(ind, 1)), kmtx(ind, 2), '-', 'Color', coll, 'LineWidth', lnw2);
plot(log10(kmtx(ind, 1)), kmtx(ind, 3), ':', 'Color', coll, 'LineWidth', lnw2);
plot(log10(kmtx(ind, 1)), kmtx(ind, 4), '-', 'Color', colr, 'LineWidth', lnw2);
plot(log10(kmtx(ind, 1)), kmtx(ind, 5), ':', 'Color', colr, 'LineWidth', lnw2);

if (tstnum <= 50)

plot(log10(kmtx(ind, 1)), kmtx(ind, 2), '.', 'Color', colb, 'Markersize', mrsz);
plot(log10(kmtx(ind, 1)), kmtx(ind, 3), '.', 'Color', colb, 'Markersize', mrsz);
plot(log10(kmtx(ind, 1)), kmtx(ind, 4), '.', 'Color', colb, 'Markersize', mrsz);
plot(log10(kmtx(ind, 1)), kmtx(ind, 5), '.', 'Color', colb, 'Markersize', mrsz);

end

xlabel('log_{10}(k_{on})  (m^3/mol)');
ylabel('Relative error  (%)');

tit = [ 'Estimation errors for K_d = ', num2str(Kd, '%.2E') ];
tit = [ tit, '  (v^* = ', num2str(max(v), '%.2E'), ', D^* = ', num2str(max(D), '%.2E'), ', l = ', num2str(l00), ', c_{eq}^* = ', num2str(max(u00), '%.2E'), ', t_{max} = ', num2str(tmax), ', x_{det} = ', num2str(xmax), ', t_{div} = ', num2str(I), ')' ];
th = title(tit);

hold off;

legend('My k_{on} error', 'My k_{off} error', 'Area k_{on} error', 'Area k_{off} error');
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(gca, 'Color', 'None');
set(gca, 'TickDir', 'out');
glg = findobj(gcf, 'Type', 'axes', 'Tag', 'legend');
set(glg, 'Color', 'None');
glg2 = findall(gcf, 'tag', 'legend');
set(glg2, 'location', 'northeastoutside', 'Color', 'None');

disp(' ');
expfig = input('Export figure? y/n = 1/0 ');
if (expfig)
	if (exist([ pwd, '\Exports']) == 0)
		system('mkdir Exports');
	end
	export_fig './Exports/estimationerrorsk' -png -transparent -m3;
end

sad = input('Export error data? y/n = 1/0 ');
if (sad)

if (exist([ pwd, '\Exports']) == 0)
	system('mkdir Exports');
end

fnm = 'estimationdatak';
xlswrite([ '.\Exports\', fnm ], {'log_10(k_on)'}, 1, 'A1');
xlswrite([ '.\Exports\', fnm ], log10(kmtx(ind, 1)), 1, 'A2');
xlswrite([ '.\Exports\', fnm ], {'R.e. k_on (mine)'}, 1, 'B1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 2), 1, 'B2');
xlswrite([ '.\Exports\', fnm ], {'R.e. k_off (mine)'}, 1, 'C1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 3), 1, 'C2');

xlswrite([ '.\Exports\', fnm ], {'R.e. k_on (area)'}, 1, 'D1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 4), 1, 'D2');
xlswrite([ '.\Exports\', fnm ], {'R.e. k_off (area)'}, 1, 'E1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 5), 1, 'E2');

end

else

clf;
hold on;

axis([ log10(k0(1)*(10^(-magfac))), log10(k0(1)*(10^magfac)), 0, 100 ]);

for i = 1:(lkp - 2)
	plot(log10(kmtx(ind, 1)), kmtx(ind, 1+i), '-',  'Color', [rand(1), rand(1), rand(1)], 'LineWidth', 2);
end

xlabel('log_{10}(k_{on})  (m^3/mol)');
ylabel('Relative error  (%)');

tit = [ 'Estimation errors for K_d = ', num2str(Kd, '%.2E') ];
tit = [ tit, '  (v^* = ', num2str(max(v), '%.2E'), ', D^* = ', num2str(max(D), '%.2E'), ', l = ', num2str(l00), ', c_{eq}^* = ', num2str(max(u00), '%.2E'), ', t_{max} = ', num2str(tmax), ', x_{det} = ', num2str(xmax), ', t_{div} = ', num2str(I), ')' ];
th = title(tit);

hold off;

legend('c_L', 's_L^1', 's_L^2', 'h_L', 'c_T', 's_T^1', 's_T^2', 'h_T', 'c_C', 's_C^1', 's_C^2', 'h_C');

set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
set(gca, 'Color', 'None');
set(gca, 'TickDir', 'out');
glg = findobj(gcf, 'Type', 'axes', 'Tag', 'legend');
set(glg, 'Color', 'None');
glg2 = findall(gcf, 'tag', 'legend');
set(glg2, 'location', 'northeastoutside', 'Color', 'None');

disp(' ');
expfig = input('Export figure? y/n = 1/0 ');
if (expfig)
	if (exist([ pwd, '\Exports']) == 0)
		system('mkdir Exports');
	end
	export_fig './Exports/estimationerrorsp' -png -transparent -m3;
end

sad = input('Export error data? y/n = 1/0 ');
if (sad)

if (exist([ pwd, '\Exports']) == 0)
	system('mkdir Exports');
end

fnm = 'estimationdatap';
xlswrite([ '.\Exports\', fnm ], {'log_10(k_on)'}, 1, 'A1');
xlswrite([ '.\Exports\', fnm ], log10(kmtx(ind, 1)), 1, 'A2');

xlswrite([ '.\Exports\', fnm ], {'R.e. c_L'}, 1, 'B1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 2), 1, 'B2');
xlswrite([ '.\Exports\', fnm ], {'R.e. s_L^1'}, 1, 'C1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 3), 1, 'C2');
xlswrite([ '.\Exports\', fnm ], {'R.e. s_L^2'}, 1, 'D1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 4), 1, 'D2');
xlswrite([ '.\Exports\', fnm ], {'R.e. h_L'}, 1, 'E1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 5), 1, 'E2');

xlswrite([ '.\Exports\', fnm ], {'R.e. c_T'}, 1, 'F1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 6), 1, 'F2');
xlswrite([ '.\Exports\', fnm ], {'R.e. s_T^1'}, 1, 'G1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 7), 1, 'G2');
xlswrite([ '.\Exports\', fnm ], {'R.e. s_T^2'}, 1, 'H1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 8), 1, 'H2');
xlswrite([ '.\Exports\', fnm ], {'R.e. h_T'}, 1, 'I1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 9), 1, 'I2');

xlswrite([ '.\Exports\', fnm ], {'R.e. c_C'}, 1, 'J1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 10), 1, 'J2');
xlswrite([ '.\Exports\', fnm ], {'R.e. s_C^1'}, 1, 'K1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 11), 1, 'K2');
xlswrite([ '.\Exports\', fnm ], {'R.e. s_C^2'}, 1, 'L1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 12), 1, 'L2');
xlswrite([ '.\Exports\', fnm ], {'R.e. h_C'}, 1, 'M1');
xlswrite([ '.\Exports\', fnm ], kmtx(ind, 13), 1, 'M2');

end

end

res = 1;