%____________________________________________________________________________________________
% Algorithm 		Flower Pollination Algorithm (FPA) by X.-S. Yang, M. Karamanoglu, X.S. He
% Source			Flower Pollination Algorithm, MathWorks.com
% Author			Xin-She Yang, May 2012
% Modified by 		József Vass at the Krylov Lab, York University
% Last revised		Aug. 26, 2016 (based on Jan. 2014 version by Yang)
% Publications		[1] X.-S. Yang, Flower pollination algorithm for global optimization
%					[2] X.-S. Yang, M. Karamanoglu, X.S. He,
%					    Multi-objective flower algorithm for optimization
%					http://arxiv.org/abs/1312.5673
%					http://arxiv.org/abs/1404.0695
% 					https://en.wikipedia.org/wiki/List_of_metaphor-based_metaheuristics#Flower_pollination_algorithm_.28Yang_2012.29
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function [ best, fmin, convd ] = OptAlgGen_FlowerPollinationLog(Fun0, xy, magfac, Tol, miter)

% Initialization:
Fun = @(x) Fun0(10.^x);
sz = size(xy);
d = sz(2) - 1; % Number of optimization variables.
%Lb = log10(xy(1, 1:d)/magfac); % Lower bounds of the search region.
%Ub = log10(xy(1, 1:d)*magfac); % Upper bounds of the search region.
Lb = log10(xy(1, 1:d));
Ub = Lb;
Lb(1) = Lb(1) - log10(magfac); % Lower bounds of the search region.
Ub(1) = Ub(1) + log10(magfac); % Upper bounds of the search region.
if (d == 2)
	Lb(2) = Lb(2) - log10(magfac);
	Ub(2) = Ub(2) + log10(magfac);
else % AG plug parameters.
	Lb(2:d) = Lb(2:d) - 2;
	Ub(2:d) = Ub(2:d) + 2;
end
n = sz(1); % Number of flowers is the number of rows in xy.
Sol = log10(xy(1:n, 1:d));
%Fitness = zeros(1, n);
%Fitness(:) = xy(:, 3);
Fitness = xy(:, d+1);
best = Sol(1, :);
fmin = xy(1, d+1);
S = Sol;
p = 0.8; % Probabibility switch.

% Recorder:
[ y, ind ] = sort(-xy(:, d+1));

% Iteration:
iter = 0; % Counter.
while ((fmin > Tol) && (iter < miter))

iter = iter + 1;

% Loop over all solutions
for i = 1:n
% Pollens are carried by insects and thus can move in large scale, large distance.
% This L should replace by Levy flights. Formula: x_i^{t+1} = x_i^t + L*(x_i^t - gbest).
if (rand > p)
	% L=rand;
	L = Levy(d);
	dS = L.*(Sol(i,:) - best);
	S(i,:) = Sol(i,:) + dS;
	S(i,:) = simplebounds(S(i,:), Lb, Ub); % Check if the simple limits/bounds are OK.
else % If not, then local pollenation of neighbor flowers.
	epsilon = rand;
	JK = randperm(n); % Find random flowers in the neighbourhood.
	% As they are random, the first two entries also random. If the flower are the same or similar species,
	% then they can be pollenated, otherwise, no action. Formula: x_i^{t+1} = x_i^t + epsilon*(x_j^t - x_k^t).
	S(i,:)=S(i,:)+epsilon*(Sol(JK(1),:)-Sol(JK(2),:));
	S(i,:)=simplebounds(S(i,:),Lb,Ub); % Check if the simple limits/bounds are OK.
end

Fnew = Fun(S(i,:)); % Evaluate new solutions.
% If fitness improves (better solutions found), update then:
if (Fnew <= Fitness(i))
	Sol(i,:) = S(i,:);
	Fitness(i) = Fnew;
end

% Update the current global best:
if (Fnew <= fmin)

best = S(i,:);
fmin = Fnew;

end % if

end % for

end % while

best = 10.^best;
convd = (fmin <= Tol);
%____________________________________________________________________________________________
% Sub-functions
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
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
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function L = Levy(d)
% Draw n Levy flight sample.
% Levy exponent and coefficient. For details, see Chapter 11 of the following book:
% Xin-She Yang, Nature-Inspired Optimization Algorithms, Elsevier, (2014).

beta = 3/2;
sigma = (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u = randn(1,d)*sigma;
v = randn(1,d);
step = u./abs(v).^(1/beta);
L = 0.01*step;