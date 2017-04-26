function [ kon, koff, kppars, frerr, convd, succi ] = Inverter_Master(Metric, LpC0, relI, kE, mth1, mth2, ptp, I, J, magfac, L, useL, acc, rerr, miter)

global ioptvec
global ioptvecmlt
global iopttime

convd = 0;
n = 0;
frerr = 10^50;
succi = 1;

while ((~convd) && (n < length(ioptvec)))

n = n + 1;
mltn = ioptvecmlt(n);
N = 0;

while ((~convd) && (N < mltn))

N = N + 1;

if (ioptvecmlt(n) == 1)
	cntstr = '';
else
	cntstr = [ num2str(N), ' ' ];
end
disp([ 'Execution ', cntstr, 'of the ', InverterName(ioptvec(n)), ' inverter for ~', num2str(iopttime(ioptvec(n))/60, '%.2f'), ' mins until ~', Min2Time(iopttime(ioptvec(n))/60), '...' ]);

[ kon0, koff0, kppars, frerr0, convd ] = Inverter(ioptvec(n), Metric, LpC0, relI, kE, mth1, mth2, ptp, I, J, magfac, L, useL, acc, rerr);

if ((n == 1) || (frerr0 < frerr))
	kon = kon0;
	koff = koff0;
	frerr = frerr0;
	if (convd)
		succi = ioptvec(n);
	end
end

end % while 2
end % while 1

%disp('Finished running the inverter(s).');