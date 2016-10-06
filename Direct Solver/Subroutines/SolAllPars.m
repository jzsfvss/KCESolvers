function g = SolAllPars(efnm, mth1, mth2, ptp)

global k % k = (k_on, k_off) vector (m^3/mol, none).
global v % v = (v_L, v_T, v_C) mobilities (m^2/V*sec).
global D % D = (D_L, D_T, D_C) diffusion coefficients (m^2/sec).
global l % l = length of the injected plug (m).
global l2 % l2 = length of the second plug (m).
global u0 % u0 = (L0eq, T0eq, C0eq) initial concentrations (mol/m^3).
global u022 % u022 = T0eq_2 initial concentration (mol/m^3).
global xmax % xmax = total capillary length (m).
global tmax % tmax = total time of run (sec).

% Importing:
T = readtable([ pwd, '\', efnm ]);
sz = size(T);
Jcol = sz(2)-3;

% Determining mesh size:
I = 0;
J = 1;
for j = 1:Jcol
	I = max(T{17, 3+j}, I);
end
if (I == 0)
	disp(' ');
	I = input('Mesh divisions in t = ');
end

% Initializing output matrix 1:
if (mth1 == 2)
	Cind = 2;
else
	Cind = 3;
end
U1 = zeros(I+1, Jcol);
U2 = zeros(I+1, Jcol*Cind);

% Generating solutions for each column of parameters:
disp(' ');
for j = 1:Jcol

% Progress indicator:
disp([ 'Generating solution for parameter column ', T.Properties.VariableNames{3+j}, '...' ]);

% Reading data:
k = T{1:2, 3+j};
v = T{3:5, 3+j};
D = T{6:8, 3+j};
l = T{9, 3+j};
l2 = T{10, 3+j};
u0 = T{11:13, 3+j};
u022 = T{14, 3+j};
xmax = T{16, 3+j};
tmax = T{15, 3+j};
if (mth1 == 2)
		v = [ v(1), v(3) ]';
		D = [ D(1), D(3) ]';
end

% Generating solution:
[ u, t, x, st, sx ] = DirSolver(mth1, mth2, ptp, I, J);
U1(:, j) = permute(u(1, :, sx+1) + u(Cind, :, sx+1), [ 2, 3, 1 ]);
for n = 1:Cind
	U2(:, (j-1)*Cind + n) = permute(u(n, :, sx+1), [ 2, 3, 1 ]);
end

end % for j

% Creating file names:
if (exist([ pwd, '\Exports']) == 0)
	system('mkdir Exports');
end
if (mth1 == 1)
	fnm = MethodName(mth2+2);
else
	fnm = MethodName(mth1-1);
end
fnm1 = [ fnm, '_', DensityName(ptp), '_', num2str(I), 'x', num2str(sx), '_L+C' ];
if (mth1 == 2)
	lbls = 'LC';
else
	lbls = 'LTC';
end
fnm2 = [ fnm, '_', DensityName(ptp), '_', num2str(I), 'x', num2str(sx), '_', lbls ];

% Exporting L+C:
header1 = T.Properties.VariableNames(4:(3+Jcol));
xlswrite([ '.\Exports\SolDS_', fnm1 ], header1);
xlswrite([ '.\Exports\SolDS_', fnm1 ], U1, 1, 'A2');

% Exporting solution components:
header2 = cell(1, Jcol*Cind);
for i = 1:Jcol
for j = 1:Cind
	header2((i-1)*Cind + j) = {[ lbls(j), '_', T.Properties.VariableNames{3+i} ]};
end
end
xlswrite([ '.\Exports\SolDS_', fnm2 ], header2);
xlswrite([ '.\Exports\SolDS_', fnm2 ], U2, 1, 'A2');

disp(' ');
disp('Solutions saved to files:');
disp([ pwd, '\Exports\SolDS_', fnm1, '.*' ]);
disp([ pwd, '\Exports\SolDS_', fnm2, '.*' ]);

g = 1;