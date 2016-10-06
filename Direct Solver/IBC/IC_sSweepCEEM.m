function ic = IC_sSweepCEEM(x, n)

global l
global u0
global u022

ic = [ u0(1)*MyHeaviside(x-l), u0(2)*Density(x/l, n) + u022*MyHeaviside(x-l), u0(3)*MyHeaviside(x-l) ];