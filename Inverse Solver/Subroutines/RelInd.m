function [ rI, rIlp, rIbr, rIrp ] = RelInd(LpC, ritol, ptp, iopt)

global tmax

% Time-derivative:
I = length(LpC)-1;
dt = tmax/I;
t = linspace(0, tmax, I+1);
dtLpC = (LpC(2:(I+1)) - LpC(1:I))/dt;

% All relevant indices:
mxLpC = max(LpC);
rI = find(LpC > ritol*mxLpC);
rI1 = min(rI);
rI2 = max(rI);
rI = rI1:rI2;
%hold on
%plot(t(rI), dtLpC(rI));
%hold off
%return

if ((ptp ~= 6) || (iopt == 4) || (iopt == 11)) % Non-AG plug, or AG with Cuckoo Search.

% Finding the bridge indices:
midind = round((rI1 + rI2)/2);
infs = zeros(1, 4);
[ lm, lmi ] = findpeaks(dtLpC(rI), 'MinPeakHeight', ritol*mxLpC);
lmil = length(lmi);
infs(1) = rI1 + lmi(1) - 1;
infs(3) = rI1 + lmi(lmil) - 1;
[ lm, lmi ] = findpeaks(-dtLpC(rI), 'MinPeakHeight', ritol*mxLpC);
lmil = length(lmi);
infs(2) = rI1 + lmi(1) - 1;
infs(4) = rI1 + lmi(lmil) - 1;
pk1 = (infs(1)+infs(2))/2;
pk2 = (infs(3)+infs(4))/2;
%fac = (rI2 - rI1)/(2*(pk1 - rI1));
fac = 4;
rIbr = [];
while (length(rIbr) < 0.4*length(rI))
	rIbr10 = round(pk1 + fac*(pk1 - rI1));
	rIbr20 = round(pk2 - fac*(rI2 - pk2));
	rIbr1 = min(rIbr10, rIbr20);
	rIbr2 = max(rIbr10, rIbr20);
	rIbr = rIbr1:rIbr2;
	fac = sqrt(fac);
end

% Setting the plug indices:
rIlp = rI1:(rIbr1-1);
rIrp = (rIbr2+1):rI2;
% Test plot:
%{
hold on
plot(t(rIlp), LpC(rIlp), 'r-');
plot(t(rIbr), LpC(rIbr), 'b-');
plot(t(rIrp), LpC(rIrp), 'r-');
hold off
%}

else % AG plug.

rIlp = [];
rIbr = [];
rIrp = [];

end