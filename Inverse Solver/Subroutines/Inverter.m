function [ kon, koff, kppars, frerr, convd, succi ] = Inverter(iopt, Metric, LpC0, relI, kE, mth1, mth2, ptp, I, J, magfac, L, useL, acc, rerr);

global iopttime
global numtdiv
global rdirsol
global optstinds

global mLpC0
global krec
global kpparsrec
global kpparsest
global l0

miter = round(iopttime(iopt)/(numtdiv(iopt)*rdirsol));
sz = size(krec);
optstinds = [ optstinds, sz(1) + 1 ];

% Method:
switch (iopt)
case 2
	oag = 1;
case 4
	oag = 2;
case 7
	oag = 4;
case 9
	oag = 3;
case 12
	oag = 5;
case 13
	oag = 6;
case 14
	oag = 7;
case 15
	oag = 8;
otherwise
	oag = 0;
end

% Solution:
if (iopt == 0) % Run all inverters.

[ kon, koff, kppars, frerr, convd, succi ] = Inverter0(Metric, LpC0, relI, kE, mth1, mth2, ptp, I, J, magfac, L, useL, acc, rerr, miter)

else % Run a specific inverter.

% Initialization:
if (useL)
	yeps = (10^(-(acc+1)))/(L*mLpC0);
	disp([ 'Threshold: relative error < ', num2str(100*yeps, '%.2E') ]);
else
	yeps = rerr;
end

% Optimization:
if (ptp ~= 6) % Standard plugs.

EF = @(k) Error(Metric, LpC0, relI, k, mth1, mth2, ptp, I, J)/mLpC0;
if (iopt == 2)
	[ k, Ek, convd ] = OptAlg_NelderMead(EF, kE, yeps, miter);
else
	[ k, Ek, convd ] = OptAlgGen_Master(oag, EF, kE, magfac, yeps, miter);
end
kon = k(1);
koff = k(2);
kppars = k;

else % AG plug.

if (mth1 ~= 2) % Standard methods.

cT = kpparsest(7);
sT = kpparsest(8:9);
EF = @(kpp) Error(Metric, LpC0, relI, KPparsTrans(kpp, cT, sT), mth1, mth2, ptp, I, J)/mLpC0;
optinds = [ 2, 3:6, 7, 11:14 ];

sz = size(kE);
lkE = sz(1);
ykE = zeros(lkE, 1);
for i = 1:lkE
	ykE(i) = EF(kE(i, optinds));
end
[ kpp, Ek, convd ] = OptAlgGen_Master(oag, EF, [ kE(:, optinds), ykE ], magfac, yeps, miter);

kppars = KPparsTrans(kpp, cT, sT);
kon = kppars(1);
koff = kppars(2);

else % MASKE.

EF = @(kpp) Error(Metric, LpC0, relI, kpp, mth1, mth2, ptp, I, J)/mLpC0;
[ kpp, Ek, convd ] = OptAlgGen_Master(oag, EF, kE, magfac, yeps, miter);
kon = kpp(1);
koff = kpp(2);
kppars = kpp;

end

end
frerr = Ek;
succi = iopt;

% Global minimum:
if (~convd) % 1
if (ptp ~= 6) % 2

[ frerr2, ind2 ] = min(krec(:, 3));
if (frerr2 < frerr) % 3
	kon = krec(ind2, 1);
	koff = krec(ind2, 2);
	frerr = krec(ind2, 3);
	convd = (frerr <= yeps);
end % 3

else % 2

sz = size(kpparsrec);
szc = sz(2);
[ frerr2, ind2 ] = min(kpparsrec(:, szc));
if (frerr2 < frerr) % 4
	kon = kpparsrec(ind2, 1);
	koff = kpparsrec(ind2, 2);
	kppars = kpparsrec(ind2, 1:(szc-1));
	frerr = kpparsrec(ind2, szc);
	convd = (frerr <= yeps);
end % 4

end % 2
end % 1

end