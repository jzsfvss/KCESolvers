function y = Error(Metric, LpC0, relI, k0, mth1, mth2, ptp, I, J)

global k
global krec
global kppars
global kpparsrec
global u0
global Lini
global Tini
global mLpC0

% Initialization:
k = [ k0(1), k0(2) ]';
kppars = k0;
u0 = ConcIni(Lini, Tini, k(2)/k(1));

% Solution generation:
LpC = SolLpC(mth1, mth2, ptp, I, J);
y = Metric(LpC0(relI), LpC(relI));

% Recording parameters:
krec = [ krec; k(1), k(2), y/mLpC0 ];
if (ptp == 6) % Asym. G. plug.
	kpparsrec = [ kpparsrec; kppars, y/mLpC0 ];
end