function [ k, v, D, l, l2, u0, u022, xmax, tmax ] = ParsInput(mth1, mth2)

kon = input('k_on (m^3/mol) = ');
koff = input('k_off (m^3/mol) = ');
k = [ kon, koff ];

vL = input('v_L (m^2/V*sec) = ');
if (mth1 ~= 2)
	vT = input('v_T (m^2/V*sec) = ');
else
	vT = 0;
end
vC = input('v_C (m^2/V*sec) = ');
v = [ vL, vT, vC ];

DL = input('D_L (m^2/sec) = ');
if (mth1 ~= 2)
	DT = input('D_T (m^2/sec) = ');
else
	DT = 0;
end
DC = input('D_C (m^2/sec) = ');
D = [ DL, DT, DC ];

l = input('Plug length (m) = ');
if ((mth1 == 1) && (mth2 == 7))
	l2 = input('Plug length 2 (m) = ');
else
	l2 = l;
end

L0eq = input('L0eq (mol/m^3) = ');
T0eq = input('T0eq (mol/m^3) = ');
if ((mth1 == 1) && (mth2 == 5))
	u022 = input('T0eq_2 (mol/m^3) = ');
else
	u022 = T0eq;
end
C0eq = input('C0eq (mol/m^3) = ');
u0 = [ L0eq, T0eq, C0eq ];

xmax = input('Detector location (m) = ');
tmax = input('Duration (sec) = ');

if (mth1 == 2)
		v = [ v(1), v(3) ];
		D = [ D(1), D(3) ];
end