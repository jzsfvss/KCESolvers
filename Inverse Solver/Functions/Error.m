function y = Error(Metric, LpC0, relI, k0, mth1, mth2, ptp, I, J, anomhand)

global k
global krec
global kppars
global kpparsrec
global u0
global Lini
global Tini
global mLpC0

%minf = Inf;
%minf = 1E100;
minf = 1E10;
%minf = 1E5;

% Initialization:
k = [ k0(1), k0(2) ]';
kppars = k0;
u0 = ConcIni(Lini, Tini, k(2)/k(1));

% Solution generation:
LpC = SolLpC(mth1, mth2, ptp, I, J);

% Error calculation:
switch (anomhand)
case 0
	y = Metric(LpC0(relI), LpC(relI));
	dorec = 1;
case 1
	if (length(LpC0) == length(LpC))
		y = Metric(LpC0(relI), LpC(relI));
		dorec = 1;
	else
		y = minf;
		dorec = 0;
	end
case 2
	if (length(LpC0) == length(LpC))
		[ lm, lmi ] = findpeaks(LpC(relI));
		if (length(lm) <= 3)
			y = Metric(LpC0(relI), LpC(relI));
			dorec = 1;
		else
			y = minf;
			dorec = 0;
		end
	else
		y = minf;
		dorec = 0;
	end
end

% Recording parameters:
if (dorec)
	krec = [ krec; k(1), k(2), y/mLpC0 ];
	if (ptp == 6) % Asym. G. plug.
		kpparsrec = [ kpparsrec; kppars, y/mLpC0 ];
	end
end