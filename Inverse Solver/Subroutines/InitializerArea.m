function [ kon, koff ] = InitializerArea(LpC)

global tmax
global Lini
global Tini

% Initialization:
I = length(LpC)-1;
dt = tmax/I;
t = linspace(0, tmax, I+1);

G = 1; % Quantum yield.
pktol = 0.01; % Peak tolerance for findpeaks.
lfac = 0.01; % Lower base threshold factor.

% Finding the peaks:
[ lm, lmi, lcl ] = MyFindPeaks(LpC, pktol);

% Determining the C peak parameters:
mC = lm(1);
tCi = lmi(1);
tC = t(tCi);
ai = min(find(LpC > lfac*mC));
a = t(ai);
dtC = tC - a;
c = tC + dtC;
ci = tCi + (tCi - ai);

% Determining the L peak parameters:
nlm = length(lm);
mL = lm(nlm);
tLi = lmi(nlm);
tL = t(tLi);
bi = max(find(LpC > lfac*mL));
b = t(bi);
dtL = b - tL;
d = tL - dtL;
di = tLi - (bi - tLi);

% Setting the intervals:
if (abs(ci - ai) > 10)
	ACi = ai:ci;
else
	ci = ai + abs(bi - di);
	ACi = ai:ci;
end
ALi = di:bi;
ADi = ci:di;

% Calculating the areas:
AC = trapz(t(ACi), LpC(ACi));
AD = trapz(t(ADi), LpC(ADi));
AL = trapz(t(ALi), LpC(ALi));

% Determining the dissipation tail areas:
hL = LpC(di);
hC = LpC(ci);
ALs = (2/3)*hL*dtL;
ACs = (2/3)*hC*dtC;

% Calculating the R1 ratio:
R1 = (AL - ALs)/(G*(AC - ACs) + ACs + ALs + AD);
% R1 = (AL - ALs)/(G*(AC - ACs) + ACs + AL + AD);
%{
tLm = tL/60;
tCm = tC/60;
ALn = (AL - ALs)/tLm;
ACn = (AC - ACs)/tCm;
ASn = (ALs + ACs + AD)/tLm;
R1 = ALn/(ALn + ACn + ASn);
%}

% Calculating the R2 ratio:
R2 = 1 + (AD + ALs + ACs)/(G*(AC - ACs));

% Calculating k:
% Kd = (Tini*(1 + R1) - Lini)/(1 + (1/R1));
Kd = 1000*(Tini*(1 + R1) - Lini)/(1 + (1/R1));
% Kd = (Tini - Lini*(1 - R1))/((1/R1) - 1);
koff = (1/tC)*log(R2);
%koff = (1/tmax)*log(R2);
% kon = (koff/Kd)/1000;
kon = koff/Kd;

% Plot:
%{
clf;
hold on;
axis([ 0, tmax, 0, 1.1*max(hC, hL) ]);
plot(t, LpC, 'k-');
l = vline(a, 'r:');
l = vline(tC, 'r:');
l = vline(c, 'r:');
l = vline(d, 'b:');
l = vline(tL, 'b:');
l = vline(b, 'b:');
hold off;
%}