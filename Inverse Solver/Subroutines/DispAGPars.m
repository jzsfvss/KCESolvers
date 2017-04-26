function dd = DispAGPars(mtht, ptp, kppars0, kppars1, kppars2)

if (ptp == 6)

N0 = length(kppars1);
kppars0 = kppars0(3:N0);
kppars1 = kppars1(3:N0);
if (~isempty(kppars2))
	kppars2 = kppars2(3:N0);
end
N0 = N0 - 2;

% Original parameters:
if (mtht == 1)

tstr = num2str(kppars0(1), '%.4E');
for n = 2:N0
	tstr = [ tstr, ', ', num2str(kppars0(n), '%.4E') ];
end
disp([ 'Original: ', tstr ]);

end

% Initial estimate:
tstr = num2str(kppars1(1), '%.4E');
for n = 2:N0
	tstr = [ tstr, ', ', num2str(kppars1(n), '%.4E') ];
end
disp([ 'Estimate: ', tstr ]);

% Solution parameters:
if (~isempty(kppars2))

tstr = num2str(kppars2(1), '%.4E');
for n = 2:N0
	tstr = [ tstr, ', ', num2str(kppars2(n), '%.4E') ];
end
disp([ 'Solution: ', tstr ]);

end

% Relative error between original and solution parameters:
if (isempty(kppars2))
	kppars2 = kppars1;
end

if (mtht == 1)

if ((kppars0(1) ~= 0) && (kppars2(1) ~= 0))
	errn = 100*abs(kppars2(1) - kppars0(1))/abs(kppars0(1));
else
	errn = 0;
end
tstr = [ num2str(errn, '%.2f'), '%' ];
for n = 2:N0
	if ((kppars0(n) ~= 0) && (kppars2(n) ~= 0))
		errn = 100*abs(kppars2(n) - kppars0(n))/abs(kppars0(n));
	else
		errn = 0;
	end
	tstr = [ tstr, ', ', num2str(errn, '%.2f'), '%' ];
end
disp([ 'R. error: ', tstr ]);

end

%{
beep on
beep
beep off
%return
%}

dd = 1;

else

disp('No asymmetric Gaussian parameters to display.');
dd = 0;

end