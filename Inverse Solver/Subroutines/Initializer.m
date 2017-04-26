function [ kon, koff, kppars ] = Initializer(LpC, rI, ptp, mth1, mth2, k0)

global v
global D
global l
%global l2
global tmax
global xmax
global Lini
global Tini
%__________________________________________________________________________________________
% Initialization
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
% Setting grids:
I = length(LpC)-1;
dt = tmax/I;
t = linspace(0, tmax, I+1);
J = 2000;
x = linspace(0, xmax, J+1);
dx = xmax/J;
u = zeros(J, 1);

% Setting kppars:
kon = 0;
koff = 0;
if (ptp == 6)
	if (mth1 == 2) % MASKE.
		nmx = 2;
	else
		nmx = 3;
	end
	kppars = zeros(1, 2 + nmx*4);
else
	kppars = [ 0, 0 ];
end

% Other settings:
pktol = 0.01;
Mean = @(x) MyMean(x, 1);
%__________________________________________________________________________________________
% Estimating k_off
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
if (isempty(k0))

% Finding the relevant indices:
% 1. Finding the relevant signal part.
rI1 = min(rI);
rI2 = max(rI);
lrI = length(rI);
trI = t(rI);
LpCrI = LpC(rI);
% 2. Finding the second derivative of the logarithm of the signal. According the propagation of a single Gaussian plug, the signal is roughly e^(-koff*t)*constant where the peak is less prominent (tail). The overlap of the tails is the bridge. So we need to find the interval where the second derivative of the logarithm is near-zero, which is the part of the bridge that is relevant for koff extraction (koff-bridge).
lgLpCrI = log(LpCrI);
dt2lgLpCrI = (lgLpCrI(3:lrI) - 2*lgLpCrI(2:(lrI-1)) + lgLpCrI(1:(lrI-2)))/(dt^2);
% 3. Finding the koff-bridge.
cl = 1;
el = -20;
while (cl)
	rI2 = find(abs(dt2lgLpCrI) < (10^el)*max(abs(dt2lgLpCrI))); % Part of the 2nd log-derivative below the 10^el threshold.
	cl = (length(rI2) < 0.2*lrI); % This means that the koff-bridge is still under 20% of the full signal length, so not long enough, thus keep increasing the exponent el. The while loop stops when the koff-bridge length is just above 20%.
	el = el+1;
end
rI21 = rI1 + min(rI2) - 1;
rI22 = rI1 + max(rI2) - 1;
rI2 = rI21:rI22; % The koff-bridge: the relevant indices for koff extraction.
lrI2 = length(rI2);
trI2 = t(rI2);

% Estimate for koff with least squares regression:
% Solving the [ t, 1 ]*[ -koff, constant ] = log(L+C) matrix equation, over the koff-bridge indices.
lhs = [ t(rI2)', 1 + 0*t(rI2)' ];
rhs = log(LpC(rI2))';
sln = lhs\rhs;
koff = abs(sln(1));

% Plot:
%{
hold on
plot(trI, LpC(rI), 'b-');
plot(trI2, exp(sln(1)*trI2 + sln(2)), 'r-'); % koff estimation region (red)
hold off
%}

else

koff = k0(2);

end

if (ptp == 6)
	kppars(2) = koff;
	%return
end
%__________________________________________________________________________________________
% Estimating kppars for L and C for asymmetric Gaussian (AG) plugs
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
if (ptp == 6) % AG plug.

% Functions:
F = cell(1, 2);
F{1} = @(t, x) normpdf(x, v(1)*t, sqrt(2*D(1)*t)); % SimNECEEM L, no Leq factor.
F{2} = @(t, x) exp(-koff*t)*normpdf(x, v(nmx)*t, sqrt(2*D(nmx)*t)); % SimNECEEM C, no Ceq factor.

% Convolution matrices:
Fmat = cell(1, 2);
for n = 1:2
	Fmat{n} = zeros(I, J);
	for i = 2:(I+1)
	for j = 2:(J+1)
		Fmat{n}(i-1, j-1) = F{n}(t(i), xmax - x(j));
	end % for 3
	end % for 2
end % for 1

% Derive the IC_L(x) vector and its Gaussian parameters:
% 1. Find the bottom 40% of the right side of L, unaffected by C (pure part).
[ lm, lmi, lcl ] = MyFindPeaks(LpC, pktol);
%nlm = length(lm);
%lmL = lm(nlm)
%lmLi = lmi(nlm);
mLpC = max(LpC);
lmrI = find(lm > 0.05*mLpC);
lm = lm(lmrI);
lmi = lmi(lmrI);
nlm = length(lm);
lmL = lm(nlm);
lmLi = lmi(nlm);
%{
if (mth1 ~= 2) % Not MASKE.
	lfac = 0.05;
else
	lfac = 0.001;
end
%}
lfac = 0.05;
rIL = find((LpC > lfac*lmL) & (LpC < 0.4*lmL) & (1:(I+1) > lmLi)); % Relevant indices.
%{
length(rIL)
hold on
plot(t, LpC, 'b-');
plot(t(rIL), LpC(rIL), 'r-');
hold off
return
%}

% 2. Inverse convolute to get the L-plug concentration IC_L, since approx. L = IC_L * F_L (at the peak).
ICL = Fmat{1}(rIL, :)\(LpC(rIL + 1)'/dx);

% 3. Extract the center and height of the L-plug.
mICL = max(ICL);
%size(ICL)
rIICL = find(ICL > 0.05*mICL);
ICL = ICL(rIICL);
xL = x(rIICL);
%size(ICL)
%plot(xL, ICL);
%return
lxL = length(xL);
[ lm, lmi, lcl ] = MyFindPeaks(ICL, pktol);
if (length(lmi) > 0)
	cLi = lmi(1);
	cL = xL(cLi); % Gaussian center.
	hL = lm(1); % Gaussian height.
	Lsucc = 1;
else
	Lsucc = 0;
end

if (Lsucc == 1) % Inverse convolution successful.
if (lcl == 1) % Peak (local max) found.

% 4. The inflection points give the standard deviations of the L-plug.
derderL = zeros(1, lxL-2);
for i = 2:(lxL-1)
	rder = (ICL(i+1) - ICL(i))/(xL(i+1) - xL(i));
	lder = (ICL(i) - ICL(i-1))/(xL(i) - xL(i-1));
	derderL(i) = (rder - lder)/(xL(i) - xL(i-1)); % 2nd derivative of the L-plug.
end
[ lm, lmi ] = min(abs(derderL(1:cLi)));
s1Li = lmi(1);
s1L = abs(xL(s1Li) - cL); % Gaussian st. dev. 1.
[ lm, lmi ] = min(abs(derderL(cLi:(lxL-2))));
xL2 = xL(cLi:(lxL-2));
s2Li = lmi(1);
s2L = abs(cL - xL2(s2Li)); % Gaussian st. dev. 2.
if (s2L == 0) s2L = s1L; end

else % No peak found, just at the interval endpoints.

s1L = cL/3;
s2L = s1L;

end % if 2

else % Inverse convolution not successful.

cL = l;
hL = lmL;
s1L = cL/3;
s2L = s1L;

end % if 1

% Derive the IC_C(x) vector and its Gaussian parameters:
% 1. Find the bottom 40% of the left side of C, unaffected by L (pure part).
[ lm, lmi, lcl ] = MyFindPeaks(LpC, pktol);
lmC = lm(1);
lmCi = lmi(1);
rIC = find((LpC > lfac*lmC) & (LpC < 0.4*lmC) & (1:(I+1) < lmCi)); % Relevant indices.

% 2. Inverse convolute to get the C-plug concentration IC_C, since C = IC_C * F_C.
ICC = Fmat{2}(rIC, :)\(LpC(rIC + 1)'/dx);

% 3. Extract the center and height of the C-plug.
mICC = max(ICC);
if (mICC ~= 0)

rIICC = find(ICC > 0.05*mICC);
ICC = ICC(rIICC);
xC = x(rIICC);
lxC = length(xC);
%plot(ICC)
[ lm, lmi, lcl ] = MyFindPeaks(ICC, pktol);
cCi = lmi(1);
cC = xC(cCi); % Gaussian center.
hC = lm(1); % Gaussian height.

else

lcl = 0;
cC = cL;
hC = hL; % Or hC = 0.

end

if (lcl == 1)

% 4. The inflection points give the standard deviations of the C-plug.
derderC = zeros(1, lxC-2);
for i = 2:(lxC-1)
	rder = (ICC(i+1) - ICC(i))/(xC(i+1) - xC(i));
	lder = (ICC(i) - ICC(i-1))/(xC(i) - xC(i-1));
	derderC(i) = (rder - lder)/(xC(i) - xC(i-1)); % 2nd derivative of the C-plug.
end
[ lm, lmi ] = min(abs(derderC(1:cCi)));
s1Ci = lmi(1);
s1C = abs(xC(s1Ci) - cC); % Gaussian st. dev. 1.
[ lm, lmi ] = min(abs(derderC(cCi:(lxC-2))));
xC2 = xC(cCi:(lxC-2));
if (~isempty(lmi))
	s2Ci = lmi(1);
	s2C = abs(cC - xC2(s2Ci)); % Gaussian st. dev. 2.
else
	s2C = s1C;
end

else

s1C = cC/3; % St. dev. is approx. 1/3 the left spread of the plug.
s2C = s1C;

end % if

kppars(3:6) = max(0, [ cL, s1L, s2L, hL ]);
kppars(2 + (nmx-1)*4 + (1:4)) = max(0, [ cC, s1C, s2C, hC ]);

%{
% Parameter Set_1:
kppars(3:6) = [ 1.5000E-02, 4.0000E-03, 4.0000E-03, 2.0942E-07 ];
kppars(2 + (nmx-1)*4 + (1:4)) = [ 1.5000E-02, 4.0000E-03, 4.0000E-03, 1.8952E-07 ];
%}

end % if
%__________________________________________________________________________________________
% Estimating k_on
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
if (ptp ~= 6) % Regular plug.

% Setting vC and DC:
if (mth1 == 2)
	vC = v(2);
	DC = D(2);
else
	vC = v(3);
	DC = D(3);
end

% Setting the mean and st. dev.:
switch ptp
case 2
	m0 = 1.1;
	s0 = 0.4;
case 3
	m0 = 1.5;
	s0 = 0.4;
%{
case 6
	if (mth1 == 2)
		m0 = kppars(7);
		s0 = kppars(9);
	else
		m0 = kppars(11);
		s0 = kppars(13);
	end
%}
end

% Setting the C solution function:
% We use heuristically the exact solution formula for Simplified NECEEM with a Gaussian plug (without the concentration factor Ceq; see below). This is fine, since NECEEM, ppKCE, and MASKE signals are all very similar to Simplified NECEEM.
Cfun = @(t, koff) l*exp(-koff*t)*normpdf(xmax, m0*l + vC*t, sqrt(((s0*l)^2) + 2*DC*t));

% Estimate for Ceq:
% 1. Find the indices of the lower 40% of the left side of the C-peak, relevant for Ceq extraction (the part of C that is almost completely unaffected by L, meaning the pure part).
[ lm, lmi ] = MyFindPeaks(LpC, pktol);
rI3 = find((LpC > 0.05*lm(1)) & (LpC < 0.4*lm(1)) & (1:(I+1) < lmi(1)));
rI31 = min(rI3);
rI32 = max(rI3);
rI3 = rI31:rI32;
% 2. We get Ceq via dividing the signal values at the relevant time points t(rI3), by the approx. C function values of Sim. NECEEM, and taking the average.
sm = 0;
cnt = 0;
for i = rI3
	Cesti = Cfun(t(i), koff);
	if (Cesti ~= 0)
		sm = sm + LpC(i)/Cesti;
		cnt = cnt + 1;
	end
end
if (cnt > 0)
	Ceq = sm/cnt;
else % If we have an extreme case where no relevant Cesti values are non-zero, then do a reasonable estimate.
	Ceq = 1000*(Lini + Tini)/2;
end

% We get the estimate for kon via the formulas relating Lini, Tini, Cini, Leq, Teq, Ceq and kon, koff:
Teq = 1000*Tini - Ceq;
kon = abs((koff/Teq)*(1000*Tini - Teq)/(Teq - 1000*(Tini - Lini)));

else % AG plug.

Leq = kppars(6)*sqrt(2*pi)*Mean(kppars(4:5)); % Height of plug x 2 pi x ave. of st. devs.
if (mth1 == 2) % MASKE.
	Ceq = kppars(10)*sqrt(2*pi)*Mean(kppars(8:9));
else % Not MASKE.
	Ceq = kppars(14)*sqrt(2*pi)*Mean(kppars(12:13));
end
%Teq = 1000*Tini - Ceq;
Teq = 1000*(Tini*(l*1000)) - Ceq; % Unclear why the l*1000 factor is necessary, but works.

kon = koff*Ceq/(Leq*Teq);
kppars(1) = kon;

if (~isempty(k0))
	kon = k0(1);
	kppars(1) = kon;
end

end

% Plot:
%{
hold on
plot(trI, LpC(rI), 'b-');
plot(t(rI3), LpC(rI3), 'g-'); % Ceq / kon estimation region (green)
plot(trI2, exp(sln(1)*trI2 + sln(2)), 'r-'); % koff estimation region (red)
hold off
%}
%__________________________________________________________________________________________
% Estimating kppars for T for asymmetric Gaussian (AG) plugs
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
if ((ptp == 6) && (mth1 ~= 2)) % Not MASKE.

% Assign the T kppars components (do this after both kon and koff => Kd is estimated):
if (((mth1 == 1) && (mth2 == 1)) || (mth1 == 3)) % (Simplified) NECEEM.

cT = (cL + cC)/2;
s1T = Mean([s1L s1C]);
s2T = Mean([s2L s2C]);

else % ppKCE.

cT = cL - 3*s1L - 3*s2L; % In ppKCE, the order of plugs is T and L, and the approx. plug length is 6 times the st. dev. Here we assume that the T-plug has the same left and right st. devs. as the L-plug.
s1T = s1L;
s2T = s2L;

end % if 2

sT = Mean([s1T s2T]);
%hT = Teq/(sqrt(2*pi)*sT);
Teq0 = 1000*Tini - Ceq;
hT = l*Teq0/(sqrt(2*pi)*sT); % Unclear why the l factor is necessary, but works.

kppars(7:10) = [ cT, s1T, s2T, hT ];

end % if 1