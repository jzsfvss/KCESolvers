%__________________________________________________________________________________________
% Algorithm         Harmony Search (HS) Algorithm by Z.W. Geem, J.H. Kim, G.V. Loganathan
% Source            Harmony Search Algorithm, MathWorks.com
% Author            Mohammad Fesanghary at Louisiana State University in Sep. 2010
% Modified by       József Vass at the Krylov Lab, York University
% Last revised      Mar. 2, 2017 (based on Jan 24, 2011 version by M. Fesanghary)
% Publications      [1] Z.W. Geem, Z. Woo, J.H. Kim, G.V. Loganathan (2001).
%                   [2] Z.W. Geem (2009), ISBN-13: 978-3642001840.   
%                   https://en.wikipedia.org/wiki/Harmony_search
%                   https://sites.google.com/site/fesangharyweb/downloads
%                   https://www.mathworks.com/matlabcentral/fileexchange/28850-harmony-search-algorithm
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function [ BestGen, BestFitness, convd ] = OptAlgGen_HarmonySearchLog(Fitness0, xy, magfac, FTol, MaxItr)
% This code has been written with Matlab 7.0.
% You can modify the simple constraint handlening method using more efficient
% methods. A good review of these methods can be found in:
% Carlos A. Coello Coello, "Theoretical and numerical constraint-handling techniques used
% with evolutionary algorithms: a survey of the state of the art".

%clc
%clear;
global MaxItr;
global FTol;
%global NVAR NG NH HMS HMCR PARmin PARmax bwmin bwmax;
global NVAR HMS HMCR PARmin PARmax bwmin bwmax;
%global HM NCHV fitness PVB BW gx;
global HM NCHV fitness PVB BW;
global BestIndex WorstIndex BestFit WorstFit BestGen currentIteration;

% Input:
Fitness = @(x) Fitness0(10.^x);
sz = size(xy);
NVAR = sz(2) - 1; % number of variables
HMS = sz(1); % harmony memory size
%MaxItr = 5000; % maximum number of iterations
%PVB = [1.0, 4; 0.6, 2; 40, 80; 20, 60]; % range of variables
%Lb = xy(1, 1:NVAR)/magfac; % Lower bounds of the search region.
%Ub = xy(1, 1:NVAR)*magfac; % Upper bounds of the search region.
Lb = log10(xy(1, 1:NVAR));
Ub = Lb;
Lb(1) = Lb(1) - log10(magfac); % Lower bounds of the search region.
Ub(1) = Ub(1) + log10(magfac); % Upper bounds of the search region.
if (NVAR == 2)
    Lb(2) = Lb(2) - log10(magfac);
    Ub(2) = Ub(2) + log10(magfac);
else % AG plug parameters.
    Lb(2:NVAR) = Lb(2:NVAR) - 2;
    Ub(2:NVAR) = Ub(2:NVAR) + 2;
end
PVB = [ Lb', Ub' ];

% Other:
HMCR = 0.9; % harmony consideration rate  0 < HMCR < 1
%NG = 6; % number of ineguality constraints
%NH = 0; % number of eguality constraints
PARmin = 0.4; % minumum pitch adjusting rate
PARmax = 0.9; % maximum pitch adjusting rate
bwmin = 0.0001; % minimum bandwidth
bwmax = 1.0; % maximum bandwidth

% Initiate matrices:
HM = zeros(HMS, NVAR);
NCHV = zeros(1, NVAR);
BestGen = zeros(1, NVAR);
fitness = zeros(1, HMS);
BW = zeros(1, NVAR);
%gx = zeros(1, NG);
% warning off MATLAB:m_warning_end_without_block

% Initialization:
%{
for i = 1:HMS
    for j = 1:NVAR
        HM(i,j) = randval(PVB(j,1), PVB(j,2));
    end
    fitness(i) = Fitness(HM(i,:));
end
%}
HM = log10(xy(:, 1:NVAR));
fitness = xy(:, NVAR + 1)';

% Iteration:
currentIteration = 0;
while ((currentIteration <= MaxItr) && (min(fitness) > FTol)) % while

PAR = (PARmax - PARmin)/(MaxItr)*currentIteration + PARmin;
coef = log(bwmin/bwmax)/MaxItr;

for pp = 1:NVAR
    BW(pp) = bwmax*exp(coef*currentIteration);
end

% Improvise a new harmony vector:
for i = 1:NVAR % for

ran = rand(1);
if (ran < HMCR) % if 1: memory consideration

index = randint(1, HMS);
NCHV(i) = HM(index, i);
pvbRan = rand(1);

if (pvbRan < PAR) % if 2: pitch adjusting

pvbRan1 = rand(1);
result = NCHV(i);

if (pvbRan1 < 0.5) % if 3

result = result + rand(1)*BW(i);

if (result < PVB(i, 2))
    NCHV(i) = result;
end

else % if 3

result = result - rand(1)*BW(i);
if (result > PVB(i, 1))
    NCHV(i) = result;
end

end % if 3

end % if 2

else % if 1

NCHV(i) = randval(PVB(i, 1), PVB(i, 2)); % random selection

end % if 1

end % for

newFitness = Fitness(NCHV);
UpdateHM(newFitness);
currentIteration = currentIteration + 1;

end % while

%BestFitness = min(fitness);
BestFitness = Fitness(BestGen);
BestGen = 10.^BestGen;
convd = (BestFitness <= FTol);
%__________________________________________________________________________________________
% Sub-functions
%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function uhm = UpdateHM(NewFit)

global NVAR MaxItr HMS;
global HM NCHV BestGen fitness;
global BestIndex WorstIndex BestFit WorstFit currentIteration;

if (currentIteration == 0) % if 1

BestFit = fitness(1);

for i = 1:HMS
    if (fitness(i) < BestFit)
        BestFit = fitness(i);
        BestIndex = i;
    end
end
            
WorstFit = fitness(1);
for i = 1:HMS
    if (fitness(i) > WorstFit)
        WorstFit = fitness(i);
        WorstIndex = i;
    end
end

end % if 1

if (NewFit < WorstFit) % if 2
            
if (NewFit < BestFit)
    HM(WorstIndex, :) = NCHV;
    BestGen = NCHV;
    fitness(WorstIndex) = NewFit;
    BestIndex = WorstIndex;
else
    HM(WorstIndex, :) = NCHV;
    fitness(WorstIndex) = NewFit;
end

WorstFit = fitness(1);
WorstIndex = 1;
for i = 1:HMS
if (fitness(i) > WorstFit)
    WorstFit = fitness(i);
    WorstIndex = i;
end
end
            
end % if 2

uhm = 1;

%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function val1 = randval(Maxv, Minv)

val1 = rand(1)*(Maxv - Minv) + Minv;

%‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
function val2 = randint(Maxv, Minv)

val2 = round(rand(1)*(Maxv - Minv) + Minv);