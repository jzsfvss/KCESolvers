function [ useL, acc, rerr, magfac, numests, ndecset, ntry, ntrymin ] = RunOpts(ropts, kon1, koff1, mth1, mth2, I, J, iopt, gI, ptp)

global Lini
global Tini
global iopttime
global rtrat
global rdirsol
global numtdiv

if (ropts == 2) % Request runtime options.

disp(' ');

if IsIn(iopt, gI)
	ro = OptSel({'Require', 'Decimal places in k', 'Relative error in L+C'});
	if (ro == 1)
		useL = 1;
	else
		useL = 0;
	end
else
	useL = 0;
end

switch (iopt)
case num2cell(gI)
	if (useL)
		acc = input('Required number of decimal places = ');
		rerr = 0;
	else
		rerr = input('Required relative error (recommended < 0.001%) = ');
		rerr = rerr/100;
		acc = 0;
	end
otherwise
	acc = 0;
	rerr = 0;
end

magfacexp = input('Orders of magnitude to search around initial estimate = ');
magfac = 10^magfacexp;

switch(iopt)
case {1, 4, 9, 12, 13, 14, 15}
	numests = input('Number of nests / flowers = ');
	%numtdiv = 2*numests;
case {2, 10}
	if (ptp ~= 6)
		numests = 3;
	else
		if (mth1 ~= 2)
			numests = 15;
		else
			numests = 11;
		end
	end	
	%numtdiv = 4;
otherwise
	numests = 1;
	%numtdiv = 1;
end
%numtdiv = [ 0, numests, 1, 2*numests, 1, 1, 1, 1, numests ];
numtdiv = [ 0, 2.25*(numests/3), 1, 2*numests, 1, 1, 1, 1, numests, numests, 2*numests, 1, 1, 1, 1 ];

if ((IsIn(iopt, gI) && (iopt < 5)) || (iopt >= 9))
	ntrymin = input('Total time for generating test points (recommended > 2 mins) = ');
	ntry = round(ntrymin*60/rdirsol);
else
	ntrymin = 0;
	ntry = 0;
end
%sits = 250;

if (IsIn(iopt, gI) && (iopt > 0))
	nmin = input('Maximum time for inverter convergence (recommended < 10 mins) = ');
	iopttime(iopt) = nmin*60;
	miter = round(nmin*60/(numtdiv(iopt)*rdirsol));
else
	miter = 0;
end

if (useL)
	ndec = acc + 1;
else
	ndec = 4;
end
ndecset = [ '%.', num2str(ndec), 'E' ];

else % Default runtime options.

%disp(' ');

useL = 0;

switch (iopt)
case num2cell(gI)
	%disp('Required relative error to the signal: 0.0001%');
	rerr = (1E-3)/100;
	acc = 0;
otherwise
	acc = 0;
	rerr = 0;
end

% Search orders of magnitude around initial estimate:
%if ((iopt == 5) || (ptp == 6))
if (iopt == 5)
	magfac = 100;
else
	magfac = 10;
end

if IsIn(iopt, gI)
	ntrymin = 2*rtrat;
	%{
	if (iopt == 5)
		ntrymin = 0.5*rtrat;
	else
		ntrymin = 2*rtrat;
	end
	%}
	ntry = round(ntrymin*60/rdirsol);
else
	ntrymin = 0;
	ntry = 0;
end
%sits = 250;

switch(iopt)
case {1, 4, 9, 12, 13, 14, 15}
	numests = 10; % Run with 10 nests / flowers.
	%numtdiv = 2*numests; % Number of direct solver evaluations per iteration in the optimization algorithm.
case {2, 10}
	%numests = 3;
	if (ptp ~= 6)
		numests = 3;
	else
		if (mth1 ~= 2)
			numests = 15;
		else
			numests = 11;
		end
	end
	%numtdiv = 4;
otherwise
	numests = 1;
	%numtdiv = 1;
end
%numtdiv = [ 0, 4, 1, 2*numests, 1, 1, 1, 1, numests ];
numtdiv = [ 0, 2.25*(numests/3), 1, 2*numests, 1, 1, 1, 1, numests, numests, 2*numests, 1, 1, 1, 1 ];

if (useL)
	ndec = acc + 1;
else
	ndec = 4;
end
ndecset = [ '%.', num2str(ndec), 'E' ];

end