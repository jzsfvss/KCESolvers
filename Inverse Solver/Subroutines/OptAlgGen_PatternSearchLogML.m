function [ xopt, yopt, convd ] = OptAlgGen_PatternSearchLogML(f0, xy, magfac, yeps, miter)

% Initialization:
sz = size(xy);
dim = sz(2) - 1;
x0 = xy(1, 1:dim);
%A = [ 1, 0; 0, 1; -1, 0; 0, -1 ];
%A = eye(dim);
%b = [ magfac*kon0; magfac*koff0; -kon0/magfac; -koff0/magfac ];
A = [];
b = [];
Aeq = [];
beq = [];
%lb = zeros(1, dim);
%ub = lb + Inf;
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
options = optimoptions('patternsearch', 'MaxFunctionEvaluations', miter, 'Display', 'off', 'StepTolerance', nz, 'MeshTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'TolBind', nz);
%options = optimoptions('patternsearch', 'ObjectiveLimit', yeps, 'MaxFunctionEvaluations', miter, 'Display', 'off', 'StepTolerance', nz, 'MeshTolerance', nz, 'FunctionTolerance', nz, 'ConstraintTolerance', nz, 'TolBind', nz);
%options = optimoptions('patternsearch', 'Display', 'off', 'FunctionTolerance', yeps, 'MaxTime', iopttime(7), 'MaxIterations', Inf, 'MaxFunctionEvaluations', Inf, 'StepTolerance', 0, 'MeshTolerance', 0);
%options = optimoptions('patternsearch', 'FunctionTolerance', yeps);

% Optimization:
[ xopt, yopt ] = patternsearch(f, x0, A, b, Aeq, beq, lb, ub, nonlcon, options);
xopt = 10.^xopt;
convd = (yopt <= yeps);