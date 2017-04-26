function opt = OptSel(opts)

nopt = max(size(opts))-1;
dl = 1;

if (nopt <= 5)
	sopt = [ opts{1, 1}, ': (1) ', opts{1, 2} ];
	for k = 3:(nopt+1)
		sopt = [ sopt, ', (', num2str(k-1), ') ', opts{1, k} ];
	end
	sopt = [ sopt, '.' ];
end

while (dl)

if (nopt <= 5)
	disp(sopt);
else
	disp([ opts{1, 1}, ':' ]);
	for k = 2:(nopt+1)
		disp([ '(', num2str(k-1), ') ', opts{1, k} ]);
	end
end	

opt = input('Selection = ');
if (sum(opt == (1:nopt)) == 0)
	disp('Improper selection! Please try again.');
	if (nopt > 5)
		disp(' ');
	end
else
	dl = 0;
end

end