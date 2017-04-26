%__________________________________________________________________________________________
% Program name		KCE Direct Solver
% Language			MatLab
% Author			József Vass <jvass@yorku.ca>
% Employer			Sergey N. Krylov <skrylov@yorku.ca>
% Institution		York University
% Version date		March 10, 2017
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%__________________________________________________________________________________________
% Global constants
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
global k % k = (k_on, k_off) vector (m^3/mol, none).
global v % v = (v_L, v_T, v_C) mobilities (m^2/V*sec).
global D % D = (D_L, D_T, D_C) diffusion coefficients (m^2/sec).
global l % l = length of the injected plug (m).
global l2 % l2 = length of the second plug (m).
global u0 % u0 = (L0eq, T0eq, C0eq) initial concentrations (mol/m^3).
global u022 % u022 = T0eq_2 initial concentration (mol/m^3).
global xmax % xmax = total capillary length (m).
global tmax % tmax = total time of run (sec).
%__________________________________________________________________________________________
% Method selection
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
disp(' ');
disp('Direct Solver for Kinetic Capillary Electrophoresis');
disp(' ');
[ mth1, mth2, ptp ] = Method(1);
if (ptp ~= 0)
	pnm = DensityName(ptp);
else
	op = Plugs();
	return
end
%__________________________________________________________________________________________
% Parameter initialization
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
disp(' ');
optprs = OptSel({ 'Parameters', 'File', 'All in file', 'Enter' });
%disp(' ');
%optprs = OptSel({ 'Parameter units', 'Standard', 'All in file', 'Enter' });

if (optprs == 2)

if (mth1 > 1)
	if (mth2 == 3)
		disp('Only the numerical solution will be generated.');
	end
	mth2 = 2;
end

% Generate and export solutions for all parameter columns in file:
g = SolAllPars('parameters.xls', mth1, mth2, ptp);

else

disp(' ');
switch (optprs)
case 1
	[ k, v, D, l, l2, u0, u022, xmax, tmax, J, I, Lini, Tini, psnm ] = ParsFile('parameters.xls', mth1, 1);
case 3
	[ k, v, D, l, l2, u0, u022, xmax, tmax ] = ParsInput(mth1, mth2);
	J = 0;
	I = 0;
end

if (I == 0)
	disp(' ');
	I = input('Mesh divisions in t = ');
end
if (J == 0)
	if (((mth1 == 2) || (mth1 == 3)) && (mth2 == 1))
		J = 1;
	else
		J = input('Mesh divisions in x = ');
	end
end
%__________________________________________________________________________________________
% KCE
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
if (mth1 == 1)

disp(' ');
disp('Numerical direct solver...')
r1 = tic;
tst = sin(2);
r2 = toc(r1);
r3 = tic;
[ u, t, x, st, sx ] = DirSolver(mth1, mth2, ptp, I, J);
r4 = toc(r3);
runtime = r4/r2;
disp([ 'Runtime = ', num2str(MyRound(r4,2)), ' sec = ', num2str(MyRound(runtime/100000,2)), 'E5 x runtime(sin(2)).' ])
st = length(t)-1;
sx = length(x)-1;

% Solution verification via relative errors:
if (sx > 5)
	[ rer, wrer, crer, lrer, lwrer ] = SolVerifierKCE(u, x, t, 0.05, 0.5);
	disp(['Ave. rel. errors = ', num2str(MyRound(rer(1),2)), '%, ', num2str(MyRound(rer(2),2)), '%, ', num2str(MyRound(rer(3),2)), '%' ]);
	disp(['Ave. rel. rel. errors = ', num2str(MyRound(wrer(1),2)), '%, ', num2str(MyRound(wrer(2),2)), '%, ', num2str(MyRound(wrer(3),2)), '%' ]);
else
	disp('Rel. errors could not be estimated.');
end

% File name:
fnm = [ psnm, '_', MethodName(mth2+2), '_', pnm, '_', num2str(I), 'x', num2str(sx) ];
disp([ 'Solved tx-mesh size = ', num2str(I), 'x', num2str(sx) ]);

end
%__________________________________________________________________________________________
% MASKE and Simplified NECEEM
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
if ((mth1 == 2) || (mth1 == 3))
if (mth2 == 1 || mth2 == 2) % Single solution u:

disp(' ');
if (mth2 == 1)
	disp('Exact direct solver...');
else
	disp('Numerical direct solver...');
end
r1 = tic;
tst = sin(2);
r2 = toc(r1);
r3 = tic;
[ u, t, x, st, sx ] = DirSolver(mth1, mth2, ptp, I, J);
r4 = toc(r3);
runtime = r4/r2;
disp([ 'Runtime = ', num2str(MyRound(r4,2)), ' sec = ', num2str(MyRound(runtime/100000,2)), 'E5 x runtime(sin(2)).' ])
st = length(t)-1;
sx = length(x)-1;

% Solution verification via relative errors:
if (mth2 == 2)
if (mth1 == 2)
	if (sx > 5)
		[ rer, wrer, crer, lrer, lwrer ] = SolVerifierMASKE(u, x, t, 0.05, 0.5);
		disp(['Ave. rel. errors = ', num2str(MyRound(rer(1),2)), '%, ', num2str(MyRound(rer(2),2)), '%' ]);
		disp(['Ave. rel. rel. errors = ', num2str(MyRound(wrer(1),2)), '%, ', num2str(MyRound(wrer(2),2)), '%' ]);
	else
		disp('Rel. errors could not be estimated.');
	end
else
	if (sx > 5)
		[ rer, wrer, crer, lrer, lwrer ] = SolVerifierKCE(u, x, t, 0.05, 0.5);
		disp(['Ave. rel. errors = ', num2str(MyRound(rer(1),2)), '%, ', num2str(MyRound(rer(2),2)), '%, ', num2str(MyRound(rer(3),2)), '%' ]);
		disp(['Ave. rel. rel. errors = ', num2str(MyRound(wrer(1),2)), '%, ', num2str(MyRound(wrer(2),2)), '%, ', num2str(MyRound(wrer(3),2)), '%' ]);
	else
		disp('Rel. errors could not be estimated.');
	end
end
end

% File name:
if (mth2 == 1)
	numanal = 'Exact';
else
	numanal = 'Numerical';
end
if (mth1 == 2)
	fnm = [ psnm, '_', 'MASKE_', numanal, '_', pnm, '_', num2str(I), 'x', num2str(sx) ];
else
	fnm = [ psnm, '_', 'SimNECEEM_', numanal, '_', pnm, '_', num2str(I), 'x', num2str(sx) ];
end
disp([ 'Solved tx-mesh size = ', num2str(I), 'x', num2str(sx) ]);

else % Both exact u1 and numerical u2 solutions:

disp(' ');
disp('Exact direct solver...')
r1 = tic;
tst = sin(2);
r2 = toc(r1);
r3 = tic;
[ u1, t1, x1, st1, sx1 ] = DirSolver(mth1, 1, ptp, I, J);
r4 = toc(r3);
runtime = r4/r2;
disp([ 'Runtime = ', num2str(MyRound(r4,2)), ' sec = ', num2str(MyRound(runtime/100000,2)), 'E5 x runtime(sin(2)).' ])
st1 = length(t1)-1;
sx1 = length(x1)-1;
% disp('Rel. errors could not be estimated.');

disp(' ');
disp('Numerical direct solver...')
r1 = tic;
tst = sin(2);
r2 = toc(r1);
r3 = tic;
[ u2, t2, x2, st2, sx2 ] = DirSolver(mth1, 2, ptp, I, J);
r4 = toc(r3);
runtime = r4/r2;
disp([ 'Runtime = ', num2str(MyRound(r4,2)), ' sec = ', num2str(MyRound(runtime/100000,2)), 'E5 x runtime(sin(2)).' ])
st2 = length(t2)-1;
sx2 = length(x2)-1;
if (sx2 > 5) % Rel. err. of a single solution.

if (mth1 == 2)
	[ rer, wrer, crer, lrer, lwrer ] = SolVerifierMASKE(u2, x2, t2, 0.05, 0.5);
	disp(['Ave. rel. errors = ', num2str(MyRound(rer(1),2)), '%, ', num2str(MyRound(rer(2),2)), '%' ]);
	disp(['Ave. rel. rel. errors = ', num2str(MyRound(wrer(1),2)), '%, ', num2str(MyRound(wrer(2),2)), '%' ]);
else
	[ rer, wrer, crer, lrer, lwrer ] = SolVerifierKCE(u2, x2, t2, 0.05, 0.5);
	disp(['Ave. rel. errors = ', num2str(MyRound(rer(1),2)), '%, ', num2str(MyRound(rer(2),2)), '%, ', num2str(MyRound(rer(3),2)), '%' ]);
	disp(['Ave. rel. rel. errors = ', num2str(MyRound(wrer(1),2)), '%, ', num2str(MyRound(wrer(2),2)), '%, ', num2str(MyRound(wrer(3),2)), '%' ]);
end

else % Rel. err. between the exact and numerical solutions.

if (mth1 ~= 2)
	Cind = 3;
	lbls = 'LTC';
else
	Cind = 2;
	lbls = 'LC';
end

disp(' ');
disp('Error between the exact and numerical solutions:');
un1 = zeros(1, st1+1);
un2 = zeros(1, st2+1);
for n = 1:Cind
	un1(:) = u1(n, :, sx1+1);
	un2(:) = u2(n, :, sx2+1);
	mun1 = sqrt(sum(un1.*un1)/length(un1)); % L2 norm.
	err = sqrt(sum((un1-un2).*(un1-un2))/length(un1)); % ||un1-un2||_2.
	%mun1 = Metric_L2(un1, 0*un1);
	%err = Metric_L2(un1, un2);
	disp([ 'Abs. err. in ', lbls(n), ' = ', num2str(err, '%.2E') ]);
	disp([ 'Rel. err. in ', lbls(n), ' = ', num2str(100*err/mun1, '%.3f'), '%' ]);
end

end
t = t2;
x = x2;
sx = sx2;

% File name:
if (mth1 == 2)
	fnm = [ psnm, '_', 'MASKE_Both_', pnm, '_', num2str(I), 'x', num2str(sx) ];
else
	fnm = [ psnm, '_', 'SimNECEEM_Both_', pnm, '_', num2str(I), 'x', num2str(sx) ];
end
% disp([ 'Solved tx-mesh size = ', num2str(I), 'x', num2str(sx) ]);

end
end
%__________________________________________________________________________________________
% Plotting
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
disp(' ');
disp('Plotting...')
plt = ((mth1 == 1) || ((mth1 == 2) && ((mth2 == 1) || (mth2 == 2))) || ((mth1 == 3) && ((mth2 == 1) || (mth2 == 2))));
% plt = 1 if single solution, 0 otherwise.
if (plt)
	plt2 = FigPlot(plt, mth1, mth2, u, u, u, t, x, t, x, t, x);
else
	plt2 = FigPlot(plt, mth1, mth2, u1, u1, u2, t1, x1, t1, x1, t2, x2);
end
%__________________________________________________________________________________________
% Exporting
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
sf = input('Export figure? y/n=1/0 ');
if (sf)
	if (exist([ pwd, '\Exports']) == 0)
		system('mkdir Exports');
	end
	export_fig tmp -png -transparent -m3;
	movefile('tmp.png', [ '.\Exports\FigDS_', fnm, '.png' ]);
	disp('Figure saved to file:');
	disp([ pwd, '\Exports\FigDS_', fnm, '.png' ]);
	disp(' ');
end

sv = input('Export solution(s)? y/n=1/0 ');
if (sv == 1)
	eopt = input('Export solution at the detector only? y/n=1/0 ');
end
if (sv)
	if (plt == 1)
		SolExporter(eopt, plt, mth1, fnm, pnm, psnm, I, sx, u, u, u);
	else
		SolExporter(eopt, plt, mth1, fnm, pnm, psnm, I, sx, u1, u1, u2);
	end
end

end

disp(' ');
disp('Finished.')
disp(' ');