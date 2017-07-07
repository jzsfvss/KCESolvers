function [ locL, kerr ] = LocLipCon(kpp0, numtst, magfac, Metric, LpC0, relI, mth1, mth2, ptp, I, J)

global kpparsest
global mLpC0

% Initialization:
if (ptp ~= 6) % Standard plugs.

EF = @(k) Error(Metric, LpC0, relI, k, mth1, mth2, ptp, I, J, 0)/mLpC0;
k0 = kpp0(1:2);
kppcon = []; % Constant part.

else % AG plug.

if (mth1 ~= 2) % Standard methods.
	cT = kpparsest(7);
	sT = kpparsest(8:9);
	EF = @(kpp) Error(Metric, LpC0, relI, KPparsTrans(kpp, cT, sT), mth1, mth2, ptp, I, J, 0)/mLpC0;
	k0 = kpp0(2);
	kppcon = kpp0([ 3:6, 7, 11:14 ]);
else % MASKE.
	EF = @(kpp) Error(Metric, LpC0, relI, kpp, mth1, mth2, ptp, I, J, 0)/mLpC0;
	k0 = kpp0(1:2);
	kppcon = kpp0(3:length(kpp0));
end	 % if 2

end % if 1

% Lipschitz constant:
locL = 0;
lk0 = length(k0);

for k = 1:numtst

% Signal error:
k1 = k0.*(magfac.^normrnd(0, 1, 1, lk0));
kpp1 = [ k1, kppcon ];
E1 = EF(kpp1);

% k error:
if (lk0 == 2)
	nrmk = max(abs(k1 - k0)./abs(k0));
else % AG st. meth.
	kpp1t = KPparsTrans(kpp1, cT, sT);
	k1t = [ kpp1t(1), k1 ];
	k0t = kpp0(1:2);
	nrmk = max(abs(k1t - k0t)./abs(k0t));
end

% Lip. const. update:
locL = max(locL, nrmk/E1);

end % for

% Error:
kpp0 = [ k0, kppcon ];
E0 = EF(kpp0);

if (lk0 == 2)
	kerr = locL*abs(k0)*E0;
else % AG st. meth.
	kpp0t = KPparsTrans(kpp0, cT, sT);
	k0t = [ kpp0t(1), k0 ];
	kerr = locL*abs(k0t)*E0;
end