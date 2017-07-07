function res = InvSolAna()
%__________________________________________________________________________________________
% Program name		KCE Inverse Solver Analysis
% Language			MatLab
% Author			József Vass <jvass@yorku.ca>
% Employer			Sergey N. Krylov <skrylov@yorku.ca>
% Institution		York University
% Version date		June 8, 2017
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%__________________________________________________________________________________________
% Global constants
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
global k % k = (k_on, k_off) vector (m^3/mol, none).
global kppars % k and additional plug optimization parameters.
global kpparsest % Initial estimate of kppars.
global krec % Matrix recording rows of (k_on, k_off, E(k_on, k_off)).
global kpparsrec % Matrix recording kppars.
global srrerr % Required relative error of the solution signal.

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
% Grand loop
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% BigPars = parameter column.
% BigI = time divisions.
% srrerr = signal relative error.
% parcoll = pararmeters collected.

parcoll = [];

for BigPars = 1:12
for srrerrpow = 0:0.5:4

srrerr = 10^(-srrerrpow);
BigI = 300;
cntcoll = [ BigPars, BigI, srrerr ];
disp([ '[ pars, I, rerr ] = [ ', num2str(cntcoll(1)), ', ', num2str(cntcoll(2)), ', ', num2str(cntcoll(3), '%.2E'), ' ]' ]);

%{
for BigPars = 1:9
for BigI = 300:50:750
%BigPars = 8;
srrerr = 1E-4;
cntcoll = [ BigPars, BigI, srrerr ];
disp([ '[ pars, I, rerr ] = [ ', num2str(cntcoll(1)), ', ', num2str(cntcoll(2)), ', ', num2str(cntcoll(3), '%.2E'), ' ]' ]);
%}
%__________________________________________________________________________________________
% Method selection
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%{
disp(' ');
disp('Inverse Solver for Kinetic Capillary Electrophoresis');

disp(' ');
%}
%invfe = OptSel({'Inversion', 'Arbitrary accuracy', 'Estimates only', 'Estimation test'});
invfe = 1;

%{
if (invfe == 3)
	res = InitializerTest();
	disp(' ');
	disp('Finished.')
	disp(' ');
	return
end
%}
%{
if (invfe == 1)
	disp(' ');
	[ mth1, mth2, ptp ] = Method(2);
else
	mth1 = 1;
	mth2 = 1;
	ptp = 6;
end	
%}
mth1 = 1;
mth2 = 1;
ptp = 6;

%{
if (invfe == 1)
	disp(' ');
	mopt = OptSel({'Metric', 'L2', 'Maximum', 'Exp L2'});
else
	mopt = 1;
end
%}
mopt = 1;
Metric = MetricFun(mopt);

%ioptvec = [ 2, 4, 3, 5, 6, 7, 8, 9 ];
if (ptp ~= 6)
	ioptvec = [ 2, 7, 4, 9, 13, 12 ]; % Currently used optimization algorithms.
else
	%ioptvec = [ 7, 13, 14, 15, 4, 2, 9, 12 ]; % AG plug.
	ioptvec = [ 13, 17, 4, 12, 9, 7, 16, 2 ]; % AG plug.
end
ioptvecmlt = [ 3, 2, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]; % Multiplicities for running each inverter above, during Inverter_Master call.
ioptvecl = length(ioptvec);
ioptvecc = cell(1, 2 + ioptvecl);
ioptvecc{1} = 'Inverter';
ioptvecc{2} = 'All';
for n = 1:ioptvecl
	ioptvecc{2 + n} = InverterName(ioptvec(n));
end
%{
if (invfe == 1)
	disp(' ');
	iopt = OptSel(ioptvecc);
else
	iopt = 2;
end
%}
iopt = 2;
ioptvec2 = [ 1, ioptvec ];
ioptnm = InverterNameS(ioptvec2(iopt));
if (iopt > 1)
	iopt = ioptvec(iopt-1);
end

krec = [];
gI = [ 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17 ]; % Good inverters.
gI2 = [ 1, 2, 4, 9, 12 ]; % Inverters that require initial test points to be generated.
%__________________________________________________________________________________________
% Target signal generation / reading from file
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%disp(' ');
%mtht = OptSel({'Test signal', 'Generate', 'File' });
mtht = 1;

if (mtht == 1) % Generated target signal.

%{
disp(' ');
%disp('Select test parameter set:');
disp(' ');
%}
[ k, v, D, l, l2, u0, u022, xmax, tmax, J, I, Lini, Tini, psnm ] = ParsFile('../Direct Solver/parameters.xls', mth1, 2, BigPars);
kon0 = k(1);
koff0 = k(2);

if (ptp == 6) % Asym. G. plug.
	%r = 0.25; % Assign 0 for exact values.
	%nmx = 3 - (mth1 == 2);
	if (mth1 == 2)
		N = [ 1, 3 ];
	else
		N = 1:3;
	end		
	kppars0 = [ k(1), k(2) ]; % Shifted-Gaussian parameters.
	for n = N % Initializing kppars for L, T, C.
		kppars0 = [ kppars0, l*1.5, l*0.4, l*0.4, l*u0(n)/(sqrt(2*pi)*(l*0.4)) ];
	end
	kppars = kppars0;
	l0 = l;
	l = 1;
	l2 = 1;
	u0(N) = 0*u0(N) + 1;
end

if (I == 0)
	%disp(' ');
	%I = input('Mesh divisions in t = ');
	I = BigI;
	% J = input('Mesh divisions in x = ');
	J = 1;
end

LpC0 = SolLpC(mth1, mth2, ptp, I, J);
ptit = [ 'Parameter set: ', psnm, ';  Original: k_{on} = ', num2str(kon0, '%.2E'), ', k_{off} = ', num2str(koff0, '%.2E'), ';  ' ];

else % Experimental target signal.

%disp('HELLO')
T = readtable([ pwd, '\parameters_exp.xls' ]);
sz = size(T);	
if (sz(2) == 4)
	j = 1;
else
	disp(' ');
	disp('Select test parameter set:');
	disp(' ');
	disp(T(1:15, [1, 2, 4:sz(2)]));
	j = input('Select parameter set column number from which signals.xls originates = ');
end

% Reading the parameters:
psnm = T.Properties.VariableNames{3+j};
Lini = T{1, 3+j};
Tini = T{2, 3+j};
u022 = 0;
l = T{12, 3+j}/1000;
l0 = l;
l2 = T{13, 3+j}/1000;
xmax = T{15, 3+j}/100;
tmax = T{14, 3+j}*60;
v = zeros(3, 1);
v(1) = xmax/(60*T{6, 3+j});
v(2) = xmax/(60*T{7, 3+j});
v(3) = xmax/(60*T{8, 3+j});
D = T{9:11, 3+j};
if (mth1 == 2)
	v = [ v(1), v(3) ]';
	D = [ D(1), D(3) ]';
end

% Reading the signal:
LpC0all = xlsread([ pwd, '\signals.xls' ])';
%LpC0 = LpC0all(:, j);
LpC0 = LpC0all(j, :);
J = 1;
ptit = [ 'Parameter set: ', psnm, ';  ' ];
%{
I = sum(LpC0 ~= 0) - 1;
LpC0 = LpC0(1:(I+1));
%}

sz = size(find(isnan(LpC0)));
if (prod(sz) == 0)
	msszall = max(size(LpC0all));
	%I = szall(1) - 1;
	I = msszall - 1;
else
	I = min(find(isnan(LpC0))) - 2;
end
LpC0 = LpC0(1:(I+1));

end

%{
% TEST PLOT:
disp('Test plot...');
if ((ptp == 6) && (mth1 > 1))
	t = linspace(0, tmax, I+1);
	plot(t, LpC0);
	return
end
%}
%__________________________________________________________________________________________
% Relevant signal parts
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
ritol = 0.01;
[ rI, rIlp, rIbr, rIrp ] = RelInd(LpC0, ritol, ptp, iopt);

if (ptp ~= 6) % Non-AG plug.

disp(' ');
%{
if (invfe == 1)
	rIopt = OptSel({'Signal portion to invert', 'All', 'Bridge'});
else
	rIopt = 1;
end
%}
rIopt = 1;

switch (rIopt)
case 1
	relI = rI;
case 2
	relI = rIbr;
end

else % AG plug.
	
%relI = rI;
relI = 1:(I+1);

end

mLpC0 = Metric(LpC0(relI), 0*LpC0(relI));

%{
if (invfe == 1) % if 1
if (iopt == 1) % if 2
	ropts = 1;
else
	%if (ptp == 6) disp(' '); end
	disp(' ');
	ropts = OptSel({'Runtime options', 'Default', 'Enter'});
end % if 2
end % if 1
%}
ropts = 1;
%__________________________________________________________________________________________
% Initial estimate
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%if (ptp == 6)
%	kppars = InitializerAGP(LpC0, rI, rIlp, rIbr, rIrp, ptp, mth1);
%end
%disp(' ');
if (ptp == 6)
	%disp('Initializing the asymmetric Gaussian plug parameters for optimization...');
%	disp('Estimating k_on, k_off and the asymmetric Gaussian plug parameters...');
	l1 = l;
	l = l0;
else
	disp('Estimating k_on, k_off...');
end
[ kon1, koff1, kppars1 ] = Initializer(LpC0, rI, ptp, mth1, mth2, []);
kpparsest = kppars1;
if (ptp == 6) l = l1; end

if (invfe == 2)

[ kon1a, koff1a ] = InitializerArea(LpC0);

%{
disp(' ');
disp('Kinetic rate constants:');
ndecset = [ '%.4E' ];
if (mtht == 1)
	disp([ 'Original: k_on = ', num2str(kon0, ndecset), ', k_off = ', num2str(koff0, ndecset) ]);
end

disp([ 'Estimate (mine): k_on = ', num2str(kon1, ndecset), ', k_off = ', num2str(koff1, ndecset) ]);
if (mtht == 1)
	% disp([ 'Relative errors in k = ', num2str(MyRound(100*abs(kon2-kon0)/kon0, 2)), '%, ', num2str(MyRound(100*abs(koff2-koff0)/koff0, 2)), '%' ]);
	disp([ 'R. error (mine): k_on: ', num2str(100*abs(kon1 - kon0)/kon0, '%.2E'), '%, k_off: ', num2str(100*abs(koff1 - koff0)/koff0, '%.2E'), '%' ]);
end

disp([ 'Estimate (area): k_on = ', num2str(kon1a, ndecset), ', k_off = ', num2str(koff1a, ndecset) ]);
if (mtht == 1)
	% disp([ 'Relative errors in k = ', num2str(MyRound(100*abs(kon2-kon0)/kon0, 2)), '%, ', num2str(MyRound(100*abs(koff2-koff0)/koff0, 2)), '%' ]);
	disp([ 'R. error (area): k_on: ', num2str(100*abs(kon1a - kon0)/kon0, '%.2E'), '%, k_off: ', num2str(100*abs(koff1a - koff0)/koff0, '%.2E'), '%' ]);
end

disp(' ');
disp('Asymmetric Gaussian plug parameters:');
dd = DispAGPars(mtht, ptp, kppars0, kppars1, []);
disp(' ');
%}

%beep on
%beep
%beep off
return

end
%__________________________________________________________________________________________
% Runtime settings
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
r1 = tic;
for i = 1:25
	tst = sin(2);
end
rsin = toc(r1)/25;

t1 = tic;
k = [ kon1, koff1 ]';
u0 = ConcIni(Lini, Tini, koff1/kon1);
for i = 1:10
	LpCest = SolLpC(mth1, mth2, ptp, I, J);
end
t2 = toc(t1)/10;
rdirsol = t2;
rtrat = rdirsol/0.2; % 0.2 sec was the rdirsol value on the computer used for programming.
iopttime = [ 0, 1, 3, 10, 10, 5, 10, 5, 10, 5, 30, 15, 20, 10, 10, 10, 20 ]*60*rtrat;
 
if (ptp == 6)
	%iopttime = [ 1, 6, 1, 12, 1, 1, 6, 1, 6, 1, 1, 6, 6 ].*iopttime;
	iopttime = 12*iopttime;
	%iopttime(13) = 0.65*iopttime(13);
	iopttime(13) = 0.05*iopttime(13);
	iopttime(17) = 0.05*iopttime(17);
end
if (mth1 == 2)
	iopttime = 2*iopttime;
end

[ useL, acc, rerr, magfac, numests, ndecset, ntry, ntrymin ] = RunOpts(ropts, kon1, koff1, mth1, mth2, I, J, iopt, gI, ptp);
%__________________________________________________________________________________________
% Initial test point generation
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
if (ptp ~= 6) disp(' '); end
if (((IsIn(iopt, gI) && (iopt < 5)) || (iopt >= 9)) && (IsIn(iopt, gI2)))
	if (ptp == 6)
		ntryfac = 4;
		ntrymin = ntryfac*ntrymin;
		ntry = ntryfac*ntry;
	end
	disp([ 'Generating test points for ~', num2str(ntrymin, '%.2f'), ' mins until ~', Min2Time(ntrymin), '...' ]);
	%[ Lmax, Lave, Lwave ] = LipConEst_old(Metric, LpC0, relI, kon1, koff1, mth1, mth2, ptp, I, J, magfac, ntry)
	if (ptp ~= 6)
		[ kvec, Evec, Lmax, Lave, Lwave ] = IniTestPts(Metric, LpC0, relI, kon1, koff1, mth1, mth2, ptp, I, J, magfac, ntry);
		%L = Lmax;
		L = Lave;
		[ Esor, ind ] = sort(Evec);
		ksor = kvec(ind, :);
		krec = [ ksor, Esor/mLpC0 ];
		kE = [ ksor(1:numests, 1:2), Esor(1:numests)/mLpC0 ];
		sz = size(krec);
		kEnum = sz(1);
	else % AG plug.
		itpg = IniTestPtsGen(Metric, LpC0, relI, kppars1, mth1, mth2, ptp, I, J, magfac, ntry);
		sz = size(kpparsrec);
		eind = sz(2);
		Evec = kpparsrec(:, eind);
		kvec = kpparsrec(:, 1:(eind-1));
		[ Esor, ind ] = sort(Evec);
		ksor = kvec(ind, :);
		kpparsrec = [ ksor, Esor ];
		kE = [ ksor(1:numests, 1:(eind-1)), Esor(1:numests) ];
		sz = size(kpparsrec);
		kEnum = sz(1);		
		L = 0;
	end
else
	if (ptp ~= 6)
		L = 0;
		kE = [ kon1, koff1, 0 ];
		kEnum = 0;
	else % AG plug.
		L = 0;
		E1 = Error(Metric, LpC0, relI, kppars1, mth1, mth2, ptp, I, J, 1)/mLpC0;
		kE = [ kppars1, E1 ];
		kEnum = 0;
	end
end

if (useL)
	disp([ 'Local Lipschitz constant estimate = ', num2str(L, '%.2E') ]);
end
%disp([ 'Lipschitz factor of the approximating signal = ', num2str(max((10^(-log10(kon0)))*abs(kon2-kon0), (10^(-log10(koff0)))*abs(koff2-koff0))/(frerr2*mLpC0), '%.2E') ]);
%{
if (ptp == 6) % Asym. Gaussian plugs.

sz = size(kE);
nrkE = sz(1);
sz = size(kppars1);
lkppars1 = sz(2);
kE = [ kE(:, 1:2), zeros(nrkE, lkppars1 - 2), kE(:, 3) ];

for i = 1:nrkE
	kE(i, 3:lkppars1) = kppars1(3:lkppars1);
end % for

end % if
%}
%__________________________________________________________________________________________
% Inversion
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
if (iopt > 1)
	%disp([ 'Execution of the ', InverterName(iopt), ' inverter for ~', num2str(iopttime(iopt)/60, '%.2f'), ' mins until ~', Min2Time(iopttime(iopt)/60), '...' ]);
end
r3 = tic;
%kppars1 = kppars; % Initial asym. Gaussian plug parameter estimates.
[ kon2, koff2, kppars2, frerr2, convd, succi ] = Inverter(iopt, Metric, LpC0, relI, kE, mth1, mth2, ptp, I, J, magfac, L, useL, acc, rerr);
r4 = toc(r3);
runtime = r4/rsin;
%__________________________________________________________________________________________
% Local Lipschitz constant estimation
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
lipdiv = 1;
lipmin = ntrymin/lipdiv;
lippts = round(ntry/lipdiv);
%disp([ 'Estimating the local Lipschitz constant for ~', num2str(lipmin, '%.2f'), ' mins until ~', Min2Time(lipmin), '...' ]);

[ locL, kerr ] = LocLipCon(kppars2, lippts, 10, Metric, LpC0, relI, mth1, mth2, ptp, I, J);
%__________________________________________________________________________________________
% Output
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%{

disp(' ');
disp('Kinetic rate constants:');
if (mtht == 1)
	disp([ 'Original: k_on = ', num2str(kon0, ndecset), ', k_off = ', num2str(koff0, ndecset) ]);
end
disp([ 'Estimate: k_on = ', num2str(kon1, ndecset), ', k_off = ', num2str(koff1, ndecset) ]);
%disp([ 'Searched: k_on_min = ', num2str(kon1/magfac, '%.2E'), ', k_on_max = ', num2str(kon1*magfac, '%.2E') ]);
disp([ 'Solution: k_on = ', num2str(kon2, ndecset), ', k_off = ', num2str(koff2, ndecset) ]);
if (mtht == 1)
	% disp([ 'Relative errors in k = ', num2str(MyRound(100*abs(kon2-kon0)/kon0, 2)), '%, ', num2str(MyRound(100*abs(koff2-koff0)/koff0, 2)), '%' ]);
	disp([ 'Relative errors in k = ', num2str(100*abs(kon2-kon0)/kon0, '%.2E'), '%, ', num2str(100*abs(koff2-koff0)/koff0, '%.2E'), '%' ]);
else
	disp('Relative errors in k could not be estimated.');
end
disp([ 'Lipschitz error estimate ', char(177), ': ', num2str(kerr(1), '%.2E'), ', ', num2str(kerr(2), '%.2E') ]);
if (kon2 - kerr(1) < 0)
	kon2s = num2str(0);
else
	kon2s = num2str(kon2 - kerr(1), ndecset);
end
if (koff2 - kerr(2) < 0)
	koff2s = num2str(0);
else
	koff2s = num2str(koff2 - kerr(2), ndecset);
end
disp([ 'Solution intervals: k_on ', char(8712), ' [ ', kon2s, ', ', num2str(kon2 + kerr(1), ndecset), ' ], k_off ', char(8712), ' [ ', koff2s, ', ', num2str(koff2 + kerr(2), ndecset), ' ]' ]);
disp([ 'Local Lipschitz constant: ', num2str(locL, '%.2E') ]);

if (ptp == 6)
	disp(' ');
	disp('Asymmetric Gaussian plug parameters:');
	dd = DispAGPars(mtht, ptp, kppars0, kppars1, kppars2);
end

disp(' ');
disp('Conclusion:');
if ((succi > 0) && convd)
	disp([ 'Successfully inverted to the desired accuracy ', num2str(100*rerr, '%.2E'), '% with the ', InverterName(succi), ' inverter.']);
elseif ((succi == 0) && convd)
	disp(['Successfully inverted to the desired accuracy ', num2str(100*rerr, '%.2E'), '%.']);
else
	disp([ 'Could not invert to the desired accuracy ', num2str(100*rerr, '%.2E'), '% within the allocated time.' ]);
end
disp([ 'Relative error of the approximating signal = ', num2str(100*frerr2, '%.2E'), '%' ]);
disp([ 'Absolute error of the approximating signal = ', num2str(100*frerr2*mLpC0, '%.2E') ]);
disp([ 'Inverter runtime = ', num2str(MyRound(r4, 2)), ' secs = ', num2str(MyRound(r4/60, 2)), ' mins = ', num2str(MyRound(runtime/100000, 2)), 'E5 x runtime(sin)' ])
%}
sz = size(krec);
%disp([ 'Error evaluations (direct solver runs) = ', num2str(kEnum), ' + ', num2str(sz(1)-kEnum), ' = ', num2str(sz(1)) ]);

%{
beep on
beep
beep off
%}

inlipint = (kon0 < kon2 + kerr(1)) && (koff0 < koff2 + kerr(2));
outcoll = [ log10(kon0), kon0, koff0, kon1, koff1, kon2, koff2, 100*abs(kon2-kon0)/kon0, 100*abs(koff2-koff0)/koff0, kerr(1), kerr(2), locL, inlipint, 100*frerr2, sz(1), convd ];
parcoll = [ parcoll; cntcoll, outcoll ];
%__________________________________________________________________________________________
% Plotting
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
%{

disp(' ');
disp('Plotting:')

% Generating the approx. signal:
k = [ kon2, koff2 ]';
kppars = kppars2;
u0 = ConcIni(Lini, Tini, koff2/kon2);
[ u, t, x, st, sx ] = DirSolver(mth1, mth2, ptp, I, J);
if (mth1 ~= 2)
	Cind = 3;
else
	Cind = 2;
end
sz = size(u);
LpC(:) = u(1, :, sz(3)) + u(Cind, :, sz(3));

% Initialization:
if (mth1 == 1)
	mnm = MethodName(mth2+2);
else
	mnm = MethodName(mth1-1);
end
%mnm = [ mnm, '_', ioptnm ];
mnm = [ mnm, '_', DensityName(ptp) ];
if (exist([ pwd, '\Exports']) == 0)
	system('mkdir Exports');
end

% Plotting and exporting signal figure:
plf = input('Plot signal figure? y/n = 1/0 ');
if (plf)

ptit = [ ptit, 'Approximation: k_{on} = ', num2str(kon2, ndecset), ', k_{off} = ', num2str(koff2, ndecset), ';  Rel. error: ', num2str(100*frerr2, '%.2E'), '%' ];
plt = InvFigPlot(mth1, mth2, LpC0, LpC, u, t, rI, relI, ptit);

saf = input('Export signal figure? y/n = 1/0 ');
if (saf)

export_fig tmp -png -transparent -m3;
fnm = [ psnm, '_', mnm, '_', ioptnm, '_Signal' ];
movefile('tmp.png', [ '.\Exports\FigIS_', fnm, '.png' ]);
disp('Figure saved to file:');
disp([ pwd, '\Exports\FigIS_', fnm, '.png' ]);

end

end

% Plotting the error surface:
if (ptp ~= 6) % 1

ple = input('Plot the error surface? y/n = 1/0 ');
if (ple) % 2

%esp = ErrorSurfPlot(150, 1, kon1, koff1);
pdg1 = 0.1;
pdg2 = 1.5; % 1-2
sz = size(krec);
nopts = sz(1) - kEnum;
lam = min(1, 50/nopts);
pdg = pdg1*(1 - lam) + pdg2*lam;
esp = ErrorSurfPlot(300, 1, pdg, kon1, koff1, kEnum);

if (esp == 0) % 3

disp('Unable to visualize the error surface due to the layout of test points.');

else
	
saf = input('Export error surface figure? y/n = 1/0 ');
if (saf) % 4

export_fig tmp -png -transparent -m3;
fnm = [ psnm, '_', mnm, '_', ioptnm, '_Error' ];
movefile('tmp.png', [ '.\Exports\FigIS_', fnm, '.png' ]);
disp('Figure saved to file:');
disp([ pwd, '\Exports\FigIS_', fnm, '.png' ]);

end % 4

end % 3

end % 2

end % 1

% Exporting data:
sad = input('Export fitted signal data? y/n = 1/0 ');
if (sad)

fnm = [ psnm, '_', mnm, '_', ioptnm, '_Data' ];
xlswrite([ '.\Exports\SolIS_', fnm ], {[ psnm, '_ExpSig' ]}, 1, 'A1');
xlswrite([ '.\Exports\SolIS_', fnm ], LpC0', 1, 'A2');
xlswrite([ '.\Exports\SolIS_', fnm ], {[ psnm, '_FitSig' ]}, 1, 'B1');
xlswrite([ '.\Exports\SolIS_', fnm ], LpC', 1, 'B2');
%xlswrite([ '.\Exports\SolIS_', fnm ], LpC0, 1, [ 'A2:A', num2str(I+2) ]);
%xlswrite([ '.\Exports\SolIS_', fnm ], LpC, 1, [ 'B2:B', num2str(I+2) ]);
%disp('Signals saved to file:');
disp('Data saved to file:');
disp([ pwd, '\Exports\SolIS_', fnm, '.*' ]);

end

disp(' ');
disp('Finished.')
disp(' ');

%}

% End of the grand loop:
end
end

% Data export:
xlswrite('out', parcoll, 1, 'A1');

% End beep:
beep on
beep
beep off

res = 1;