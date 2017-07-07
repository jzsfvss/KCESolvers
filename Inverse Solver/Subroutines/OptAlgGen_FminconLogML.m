function [ xopt, yopt, convd ] = OptAlgGen_FminconLogML(fmcopt, f0, xy, magfac, yeps, miter)

dem = 0; % Display optimization exit message?

% Initialization:
sz = size(xy);
dim = sz(2) - 1;
x0 = xy(1, 1:dim);
A = [];
b = [];
Aeq = [];
beq = [];
lb = x0;
ub = x0;
lb(1) = lb(1)/magfac;
ub(1) = ub(1)*magfac;
if (dim == 2)
	lb(2) = lb(2)/magfac;
	ub(2) = ub(2)*magfac;
else % AG plug parameter ranges.
	lb(2:dim) = lb(2:dim)/100;
	ub(2:dim) = ub(2:dim)*100;
end
nonlcon = [];

% Log space:
f = @(x) f0(10.^x);
x0 = log10(x0);
lb = log10(lb);
ub = log10(ub);

% Options:
nz = 1E-500; % Very tiny near-zero value, to prevent premature termination.
mi = Inf;
if (dem)
	disp(' ');
	disp('Optimization exit message:');
	demopt = 'final-detailed';
else
	demopt = 'off';
end
switch (fmcopt)
case {1, 4}
	if (fmcopt == 1)
		hsnm = 'bfgs'; % Calculates the Hessian by a dense quasi-Newton approximation.
	else
		hsnm = 'lbfgs'; % Calculates the Hessian by a limited-memory, large-scale quasi-Newton approximation. The default memory, 10 iterations, is used.
	end
	options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'HessianApproximation', hsnm, 'ObjectiveLimit', yeps, 'MaxFunctionEvaluations', miter, 'MaxIterations', mi, 'StepTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'OptimalityTolerance', nz, 'Display', demopt);
	%options = optimoptions('fmincon', 'Algorithm', 'interior-point', 'MaxFunctionEvaluations', miter, 'StepTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'Display', 'final-detailed');
	%options = optimoptions('fmincon', 'MaxFunctionEvaluations', miter, 'StepTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'Display', 'off');
	%options = optimoptions('patternsearch', 'MaxFunctionEvaluations', miter, 'StepTolerance', nz, 'MeshTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'TolBind', nz, 'Display', 'off');
case 2
	options = optimoptions('fmincon', 'Algorithm', 'sqp', 'ObjectiveLimit', yeps, 'MaxFunctionEvaluations', miter, 'MaxIterations', mi, 'StepTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'Display', 'off');
	%options = optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxFunctionEvaluations', miter, 'StepTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'Display', 'off');
case 3
	options = optimoptions('fmincon', 'Algorithm', 'active-set', 'ObjectiveLimit', yeps, 'MaxFunctionEvaluations', miter, 'MaxIterations', mi, 'StepTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'Display', 'off');
	%options = optimoptions('fmincon', 'Algorithm', 'active-set', 'MaxFunctionEvaluations', miter, 'StepTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'Display', 'off');
end

% Optimization:
[ xopt, yopt ] = fmincon(f, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
xopt = 10.^xopt;
convd = (yopt <= yeps);