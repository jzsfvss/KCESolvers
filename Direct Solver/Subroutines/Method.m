function [ mth1, mth2, ptp ] = Method(di)

if (di == 1) % Direct Solver

mth1 = OptSel({ 'Method', 'Standard methods', 'MASKE', 'Simplified NECEEM', 'Plugs' });

if (mth1 == 4)

ptp = 0;
mth1 = 0;
mth2 = 0;

else

disp(' ');
switch (mth1)
case 1
	opts = { 'Sub-method', MethodName(3), MethodName(4), MethodName(5), MethodName(6), MethodName(7), MethodName(8), MethodName(9) };
case 2
	opts = { 'Sub-method', 'Exact', 'Numerical', 'Both' };
case 3
	opts = { 'Sub-method', 'Exact', 'Numerical', 'Both' };
end
mth2 = OptSel(opts);

disp(' ');
if ((mth1 == 3) && (mth2 ~= 2))
	opts = { 'Plug type', DensityName(1), DensityName(2), DensityName(3) };
else
	opts = { 'Plug type', DensityName(1), DensityName(2), DensityName(3), DensityName(4), DensityName(5) };
end
ptp = OptSel(opts);

end

else % Inverse Solver (di == 2)

mth1 = OptSel({ 'Method', 'Standard methods', 'MASKE', 'Simplified NECEEM' });

if (mth1 == 1)
	disp(' ');
	optvec = [ 1, 7 ];
	opts = { 'Sub-method', MethodName(3), MethodName(9) };
	mth20 = OptSel(opts);
	mth2 = optvec(mth20);
else
	mth2 = 2;
end

disp(' ');
opts = { 'Plug type', DensityName(2), DensityName(3) };
ptp = OptSel(opts) + 1;

end