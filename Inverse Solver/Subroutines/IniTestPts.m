function [ kall, Eall, Lmax, Lave, Lwave ] = IniTestPts(Metric, LpC0, relI, kon0, koff0, mth1, mth2, ptp, I, J, magfac, ntry)

global k
global u0
global Lini
global Tini
global mLpC0

%kon = [ 0, 0 ];
%koff = [ 0, 0 ];
n = 0;
Lmax = 0;
sL = 0;
nL = 0;
swL = 0;
nwL = 0;
kc = cell(2);
rn = cell(2);
LpC = cell(2);
dis0 = [ 0, 0 ];
sq2 = sqrt(2);

kall = [ kon0, koff0 ];
k = kall';
u0 = ConcIni(Lini, Tini, k(2)/k(1));
LpC1 = SolLpC(mth1, mth2, ptp, I, J);
Eall = Error2(Metric, LpC1, LpC0, relI);

% Orders of magnitude:
omon = round(log10(kon0));
omoff = round(log10(koff0));

while (n <= ntry-2)

rn{1} = [ normrnd(0, 1), normrnd(0, 1) ]*log10(magfac);
rn{2} = [ normrnd(0, 1), normrnd(0, 1) ]*log10(magfac);
for i = 1:2
	% kc{i} = ([ kon0, koff0 ].*[ magfac^normrnd(0, 1), magfac^normrnd(0, 1) ])';
	lkon = log10(kon0) + (rn{i}(1) - 0.25*rn{i}(2))/sq2;
	lkoff = log10(koff0) + (rn{i}(1) + 0.25*rn{i}(2))/sq2;
	kc{i} = [ 10^lkon, 10^lkoff ]';
	k = kc{i};
	u0 = ConcIni(Lini, Tini, kc{i}(2)/kc{i}(1));
	LpC{i} = SolLpC(mth1, mth2, ptp, I, J);
end

lhs = max((10^(-omon))*abs(kc{1}(1) - kc{2}(1)), (10^(-omoff))*abs(kc{1}(2) - kc{2}(2))); % Metric: d_inf(k1, k2).
rhs = Error2(Metric, LpC{1}, LpC{2}, relI);
dis0(1) = Error2(Metric, LpC{1}, LpC0, relI);
dis0(2) = Error2(Metric, LpC{2}, LpC0, relI);
kall = [ kall; kc{1}(1), kc{1}(2); kc{2}(1), kc{2}(2) ];
Eall = [ Eall; dis0(1); dis0(2) ];

if (max(dis0) ~= 0)
	wei = 1/max(dis0);
end

if ((max(dis0) ~= 0) && (rhs ~= 0) && (isnan(lhs) ~= 1) && (isnan(rhs) ~= 1))
	Lmax = max(Lmax, lhs/rhs);
	sL = sL + (lhs/rhs);
	nL = nL + 1;
	swL = swL + wei*(lhs/rhs);
	nwL = nwL + wei;
end

n = n+2;

end

if (nL ~= 0)
	Lave = sL/nL;
else
	Lave = Lmax;
end
if (nwL ~= 0)
	Lwave = swL/nwL;
else
	Lwave = Lmax;
end