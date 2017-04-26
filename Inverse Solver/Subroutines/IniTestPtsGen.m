function itpg = IniTestPtsGen(Metric, LpC0, relI, kppars1, mth1, mth2, ptp, I, J, magfac, ntry)

global mLpC0

% Initialization:
n = 0;
kpl = length(kppars1);
kon1 = kppars1(1);
koff1 = kppars1(2);
EF = @(k) Error(Metric, LpC0, relI, k, mth1, mth2, ptp, I, J)/mLpC0;
ef = EF(kppars1);

for n = 1:ntry

rn = normrnd(0, 1, 1, kpl); % Vector of kpl normally distributed random numbers.
lkon = log10(kon1) + (rn(1) - 0.25*rn(2))*log10(magfac)/sqrt(2);
lkoff = log10(koff1) + (rn(1) + 0.25*rn(2))*log10(magfac)/sqrt(2);

kc = [ 10^lkon, 10^lkoff, kppars1(3:kpl).*(magfac.^rn(3:kpl)) ]; % Current kppars.
ef = EF(kc); % The evaluation of the Error function records kppars and the error value into kpparsrec.

end % for

itpg = 1; % Output.