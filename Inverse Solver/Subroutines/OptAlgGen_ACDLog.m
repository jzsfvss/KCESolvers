%__________________________________________________________________________________________
% Algorithm           Adaptive Coordinate Descent
% Source              Fast Adaptive Coordinate Descent for Non-linear Optimization, MathWorks.com
% Authors             Ilya Loshchilov, Marc Schoenauer, Michele Sebag (2012)
% Modified by         József Vass at the Krylov Lab, York University
% Last revised        May 19, 2017 (based on Jun. 7, 2012 version by Ilya Loshchilov)
% Publications        [1] I. Loshchilov, M. Schoenauer, M. Sebag, Adaptive Coordinate Descent
%                     http://hal.inria.fr/docs/00/58/75/34/PDF/AdaptiveCoordinateDescent.pdf
% Notes               This source code includes the Adaptive Encoding procedure by N. Hansen, 2008.
% License             BSD license.
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function [ xmean, bestFit, convd ] = OptAlgGen_ACDLog(fobj0, xy, magfac, Tol, miter)
% strfitnessfct = target function
% N = dimension
% X_a = lower bound on all coordinates
% X_b = upper bound on all coordinates
% MAX_EVAL = max. no. of tar. fun. evaluations
% stopfitness = stopping threshold
% howOftenUpdateRotation = 1 (at each iteration, so quadratic time complexity of the algorithm, but need less function evaluations to reach the optimum)

% Initialization:
strfitnessfct = @(x) fobj0(10.^x);
sz = size(xy);
N = sz(2) - 1;
MAX_EVAL = miter;
stopfitness = Tol;
howOftenUpdateRotation = 1;
k_succ = 2.0;       
k_unsucc = 0.5;
c1 = 0.5/N;
cmu = 0.5/N;
howOftenUpdateRotation = floor(howOftenUpdateRotation); % Integer >= 1.
TolHistFun = 1e-25;
%nIterHist = 10 + floor(20*N^1.5);
arrFstop = 1e+30;
posArrF = 0;

% Set bounds:
Lb = log10(xy(1, 1:N));
Ub = Lb;
Lb(1) = Lb(1) - log10(magfac); % Lower bounds of the search region.
Ub(1) = Ub(1) + log10(magfac); % Upper bounds of the search region.
if (N == 2)
  Lb(2) = Lb(2) - log10(magfac);
  Ub(2) = Ub(2) + log10(magfac);
else % AG plug parameters.
  Lb(2:N) = Lb(2:N) - 2;
  Ub(2:N) = Ub(2:N) + 2;
end

% Initialize mean and sigma:
for i = 1:N
    XminPrime(i, 1) = Lb(i);
    XmaxPrime(i, 1) = Ub(i);
    xmean(i, 1) = XminPrime(i, 1) + rand()*(XmaxPrime(i, 1) - XminPrime(i, 1));
    sigma(i, 1) = (XmaxPrime(i, 1) - XminPrime(i, 1))/4;
end

bestFit = 1e+30;
neval = 0;

B = eye(N, N);
if 1 < 2  % Generate if desired an arbitrary orthogonal initial transformation. (N. Hansen)
	for i = 1:N
        v = randn(N, 1);
        for j = 1:i-1
          v = v - (v'*B(:, j))*B(:, j);
        end
        B(:, i) = v/norm(v);  % might fail with a very small probability
 	end
end

iter = 0;
firstAE = 1;
xperm = xpermutation(N);
% xperm = 1:N; % Non-random.
qix = 0;
somebetter = 0;

allx = zeros(N, 2*N);
allf = zeros(1, 2*N);
xmeanold = xmean;
%arrFstop = 1e+10*ones(nIterHist);
firstsav = 1;

% Generation loop:
while ((neval < MAX_EVAL) && (bestFit > stopfitness))

iter = iter + 1;
qix = qix + 1;
if (qix > N) qix = 1; end
ix = xperm(qix); % Will perform the search along this ix'th principal component.

% Sample two candidate solutions:
dx = sigma(ix, 1)*B(:, ix); % Shift along ix'th principal component, the computational complexity is linear.
x1 = xmean - dx; % First point to test along ix'th principal component.
x2 = xmean + dx; % Second point to test is symmetric to the first one on the same ix'principal component.

% Compute fitness:
Fit1 = feval(strfitnessfct, x1');
neval = neval + 1;
Fit2 = feval(strfitnessfct, x2');   
neval = neval + 1;

% Check stopping criteria:
if (0)
  posArrF = posArrF + 1;   
  if (posArrF == nIterHist) posArrF = 1; end
  arrFstop(posArrF) = Fit1;
  posArrF = posArrF + 1;
  if (posArrF == nIterHist) posArrF = 1; end
  arrFstop(posArrF) = Fit2;
  if (neval > nIterHist)
  if ((max(arrFstop) - min(arrFstop)) < TolHistFun)
    % disp('restart');
    return
  end
  end
end

% Which is the next mean point?  
lsucc = 0;
if (Fit1 < bestFit)
  bestFit = Fit1;
  xmean = x1;
  lsucc = 1;
end
if (Fit2 < bestFit)
  bestFit = Fit2;
  xmean = x2;
  lsucc = 1;
end

% Adapt step-size sigma depending on the success/unsuccess of the previous search:
if (lsucc == 1) % Increase the step-size.
  sigma(ix, 1) = sigma(ix, 1)*k_succ; % Default k_succ = 2.0.
  somebetter = 1;
else % Decrease the step-size.
  sigma(ix, 1) = sigma(ix, 1)*k_unsucc; % Default k_unsucc = 0.5.
end

% Update archive:
posx1 = ix*2 - 1;
posx2 = ix*2;
allx(:, posx1) = x1(:, 1);
allf(1, posx1) = Fit1;
allx(:, posx2) = x2(:, 1);
allf(1, posx2) = Fit2;

% Update the encoding:
update = 0;
step = floor(N/1);
if (qix == step) update = 1; end % We update our rotation matrix B every N = dimension iterations.
    
if (update && (somebetter == 1)) && 1 && (rand() < 1.1) 
  somebetter = 0;          
  [ allfsort, arindex ] = sort(allf, 2, 'ascend');
  popsz = step;
  allxbest = allx(:, arindex(1:popsz));
  if (firstAE == 1)
    ae = ACD_AEupdateFAST([], allxbest, c1, cmu,howOftenUpdateRotation); % Initialize encoding.
    ae.B = B;
    ae.Bo = ae.B; % Assuming the initial B is orthogonal.
    ae.invB = ae.B'; % Assuming the initial B is orthogonal.
    firstAE = 0;
  else 
    ae = ACD_AEupdateFAST(ae, allxbest, c1, cmu,howOftenUpdateRotation); % Adapt encoding. 
  end
  B = ae.B;
end

%{
if (rem(neval, 1000) == 0)
  % disp([ num2str(neval) ' ' num2str(bestFit-fgeneric('ftarget')) ]);
  disp([ num2str(neval) ' ' num2str(bestFit) ]);
end
%}

end % while

%outplot(out);
%outplot(inc_out);
%disp([ num2str(neval) ' ' num2str(bestFit) ]);

xmean = (10.^xmean)';
convd = (bestFit <= Tol);
%__________________________________________________________________________________________
% Sub-functions
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function arr = xpermutation(N)
% Permutation of [1,...,N].

arr1 = 1:N;
arr = arr1;
sz = N;

for i = 1:N
  pos = ceil(rand()*sz);
  arr1(i) = arr(pos);
  arr(pos) = arr(sz);
  arr(sz) = -1;
  sz = sz - 1;
end

arr = arr1;
%__________________________________________________________________________________________
% Subroutines
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function ae = ACD_AEupdateFAST(ae, pop, c1, cmu, howOftenUpdateRotation)

% This code is based on the original Adaptive Encoding procedure (and code) 
% proposed by N. Hansen, see PPSN X paper for details.

% ---------------------------------------------------------------
% Adaptive encoding. To be used under the terms of the BSD license
% Author : Nikolaus Hansen, 2008.  
% e-mail: nikolaus.hansen AT inria.fr 
% URL:http://www.lri.fr/~hansen
% REFERENCE: Hansen, N. (2008). Adaptive Encoding: How to Render
%    Search Coordinate System Invariant. In Rudolph et al. (eds.)
%    Parallel Problem Solving from Nature, PPSN X,
%    Proceedings, Springer. http://hal.inria.fr/inria-00287351/en/
% ---------------------------------------------------------------

if (nargin < 2 || isempty(pop))
  error('need two arguments, first can be empty');
end
N = size(pop, 1);

% Initialize object:
if isempty(ae) 

% Parameter setting:
ae.N = N; 
ae.mu = size(pop, 2);
ae.weights = ones(ae.mu, 1)/ae.mu;
ae.mucov = ae.mu; % For computing c1 and cmu.

if (11 < 3) % Non-uniform weights, assumes a correct ordering of input arguments.
  ae.weights = log(ae.mu+1)-log(1:ae.mu)'; 
  ae.weights = ae.weights/sum(ae.weights); 
end
ae.alpha_p = 1; 
ae.c1 = c1;
ae.cmu = cmu;
ae.cc = 1/sqrt(N); 

% Initialization:
ae.pc = zeros(N, 1);
ae.pcmu = zeros(N, 1);
ae.xmean = pop*ae.weights;
ae.C = eye(N); 
ae.Cold = ae.C;
ae.diagD = ones(N, 1);
ae.ps = 0;
ae.iter = 1;
return

end

% Begin adaptive encoding procedure:
ae.iter = ae.iter + 1;
ae.xold = ae.xmean; 
ae.xmean = pop*ae.weights;
% ae.xmean = pop(:,1);
updatePath = 1;

if (updatePath == 1) % Adapt the encoding.

dpath = ae.xmean-ae.xold;

if (sum((ae.invB*dpath).^2) == 0)
  z = 0;
else
  alpha0 = sqrt(N) / sqrt(sum((ae.invB*dpath).^2));
  z = alpha0 * dpath;
end

ae.pc = (1-ae.cc)*ae.pc + sqrt(ae.cc*(2-ae.cc))*z;
S = ae.pc*ae.pc';
ae.C = (1-ae.c1) * ae.C + ae.c1 * S;

end

if (((rem(ae.iter,howOftenUpdateRotation) == 0) || (ae.iter <= 2)) && 1)

if ((sum(isnan(ae.C(:)))) || sum(isinf(ae.C(:))))
  ae.C = ae.Cold;
end
ae.C = triu(ae.C) + triu(ae.C,1)';

% tic
[ae.Bo, EV] = eig(ae.C);
EV = diag(EV);
% toc/(2*N*N)
if (1) % limit condition of C to 1e14 + 1

cond = 1e14;
if min(EV) <= 0
  EV(EV < 0) = 0;
  tmp = max(EV)/cond;
  ae.C = ae.C + tmp*eye(N,N); EV = EV + tmp*ones(N,1); 
end
if (max(EV) > cond*min(EV))
  tmp = max(EV)/cond - min(EV);
  ae.C = ae.C + tmp*eye(N, N);
  EV = EV + tmp*ones(N, 1); 
end

end
  
ae.diagD = sqrt(EV); 
if (min(EV) <= 0)
  return
end

ae.B = ae.Bo * diag(ae.diagD);
ae.invB = diag(1./ae.diagD) * ae.Bo';
ae.Cold = ae.C;

end