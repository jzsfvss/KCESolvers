function [ k, v, D, l, l2, ueq, ueq22, xmax, tmax, J, I, Lini, Tini, psnm ] = ParsFile(fnm, mth1, optdi)

% Importing:
T = readtable([ pwd, '\', fnm ]);
sz = size(T);

% Displaying:
%disp(T(1:18, [1, 4:sz(2)]));
disp(T(20:41, [1, 2, 4:sz(2)]));
j = input('Select parameter set column number = ');

% Reading data:
psnm = T.Properties.VariableNames{3+j};
k = T{1:2, 3+j};
v = T{3:5, 3+j};
D = T{6:8, 3+j};
l = T{9, 3+j};
l2 = T{10, 3+j};
ueq = T{11:13, 3+j};
ueq22 = T{14, 3+j};
xmax = T{16, 3+j};
tmax = T{15, 3+j};
J = T{18, 3+j};
I = T{17, 3+j};
if (mth1 == 2)
		v = [ v(1), v(3) ]';
		D = [ D(1), D(3) ]';
end
if (optdi == 2)
	Lini = T{23, 3+j};
	Tini = T{24, 3+j};
else
	Lini = 0;
	Tini = 0;
end