%__________________________________________________________________________________________
% Algorithm 		Cuckoo Search (CS) Algorithm by Xin-She Yang and Suash Deb
% Source			Cuckoo Search (CS) Algorithm, MathWorks.com
% Author			Xin-She Yang at Cambridge University from Nov. 2008 to Jun. 2009
% Modified by 		József Vass at the Krylov Lab, York University
% Last revised		Aug. 15, 2016 (based on Dec. 2009 version by Yang)
% Publications		[1] X.-S. Yang, S. Deb, Cuckoo search via Levy flights
%					[2] X.-S. Yang, S. Deb, Engineering optimization by cuckoo search
%					[3] X.-S. Yang, Nature-Inspired Metaheuristic Algoirthms
%					http://arxiv.org/PS_cache/arxiv/pdf/1003/1003.1594v1.pdf					
%					http://arxiv.org/PS_cache/arxiv/pdf/1005/1005.2908v2.pdf 					
% 					https://en.wikipedia.org/wiki/Cuckoo_search
% Notes 			This version is more efficient than the one in [3].
%					The algorithm is heuristic and non-deterministic.
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function [ bestnest, fmin, convd ] = OptAlgGen_CuckooSearchLog(fobj0, xy, magfac, Tol, miter)
% Optimization of fobj in log-space.

% Initialization:
fobj = @(x) fobj0(10.^x);
sz = size(xy);
nd = sz(2) - 1; % Number of optimization variables.
Lb = log10(xy(1, 1:nd));
Ub = Lb;
Lb(1) = Lb(1) - log10(magfac); % Lower bounds of the search region.
Ub(1) = Ub(1) + log10(magfac); % Upper bounds of the search region.
if (nd == 2)
	Lb(2) = Lb(2) - log10(magfac);
	Ub(2) = Ub(2) + log10(magfac);
else % AG plug parameters.
	Lb(2:nd) = Lb(2:nd) - 2;
	Ub(2:nd) = Ub(2:nd) + 2;
end

n = sz(1); % Number of nests is the number of rows in xy.
nest = log10(xy(1:n, 1:nd));
pa = 0.25; % Discovery rate of alien eggs/solutions.

% Recorder:
[ y, ind ] = sort(-xy(:, nd + 1));

% Calculating the fitness of the initial solutions:
fitness = (10^10)*ones(n, 1); % Initializing fitness with a large value.
[ fmin, bestnest, nest, fitness ] = get_best_nest(fobj, nest, nest, fitness); % Current best.

% Iteration:
iter = 0; % Counter 1.
N_iter = 0; % Counter 2.
while ((fmin > Tol) && (iter < miter))

% Generate new solutions, but keep the current best:
new_nest = get_cuckoos(nest, bestnest, Lb, Ub);   
[ fnew, best, nest, fitness ] = get_best_nest(fobj, nest, new_nest, fitness);

% Update the counter:
N_iter = N_iter + n;

% Discovery and randomization:
new_nest = empty_nests(nest, Lb, Ub, pa);

% Evaluate this set of solutions:
[fnew, best, nest, fitness] = get_best_nest(fobj, nest, new_nest, fitness);

% Update the counter again:
N_iter = N_iter + n;

% Find the best objective so far:
if (fnew < fmin)

fmin = fnew;
bestnest = best;

end % if

iter = iter + 1;

end % while

bestnest = 10.^bestnest;
convd = (fmin <= Tol);
%__________________________________________________________________________________________
% Sub-functions
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function nest = get_cuckoos(nest, best, Lb, Ub)
% Get cuckoos by ramdom walk.

% Levy flights (Levy exponent and coefficient; see eqn (2.21) on pg 16 [3]):
n = size(nest, 1);
beta = 3/2;
sigma = (gamma(1 + beta)*sin(pi*beta/2)/(gamma((1 + beta)/2)*beta*2^((beta - 1)/2)))^(1/beta);

for j = 1:n

s = nest(j, :); % A simple way of implementing Levy flights (Mantegna's algorithm); for standard random walks use s = 1.
u = randn(size(s))*sigma;
v = randn(size(s));
step = u./abs(v).^(1/beta);
stepsize = 0.01*step.*(s - best);
% (s - best) means when the solution is the best solution, it remains unchanged.
% The factor 0.01 comes from the fact that L/100 should be the typical step size of walks/flights
% where L is the typical lenghtscale; otherwise, Levy flights may become too aggresive/efficient, 
% which makes new solutions (even) jump outside of the design domain (and thus wasting evaluations).

% Now the actual random walks or flights:
s = s + stepsize.*randn(size(s));

% Apply simple bounds / limits:
nest(j, :) = simplebounds(s, Lb, Ub);

end % for
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function [ fmin, best, nest, fitness ] = get_best_nest(fobj, nest, newnest, fitness)
% Find the current best nest.

% Evaluating all new solutions:
for j = 1:size(nest, 1)
	fnew = fobj(newnest(j, :));
    if (fnew <= fitness(j))
		fitness(j) = fnew;
		nest(j, :) = newnest(j, :);
    end
end

% Find the current best:
[ fmin, K ] = min(fitness);
best = nest(K, :);
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function new_nest = empty_nests(nest, Lb, Ub, pa)
% Replace some nests by constructing new solutions / nests.

n = size(nest, 1); % A fraction of worse nests are discovered with a probability pa.
K = (rand(size(nest)) > pa); % Discovered or not - a status vector.

% In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
% this cuckoo's egg is less likely to be discovered, thus the fitness should 
% be related to the difference in solutions. Therefore, it is a good idea 
% to do a random walk in a biased way with some random step sizes.  

% New solution by biased/selective random walks:
stepsize = rand*(nest(randperm(n), :) - nest(randperm(n), :));
new_nest = nest + stepsize.*K;
for j = 1:size(new_nest, 1)
	s = new_nest(j, :);
	new_nest(j, :) = simplebounds(s, Lb, Ub);  
end
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function s = simplebounds(s, Lb, Ub)
% Application of simple constraints.

% Apply the lower bound:
ns_tmp = s;
I = (ns_tmp < Lb);
ns_tmp(I) = Lb(I);

% Apply the upper bounds:
J = (ns_tmp > Ub);
ns_tmp(J) = Ub(J);

% Update this new move:
s = ns_tmp;