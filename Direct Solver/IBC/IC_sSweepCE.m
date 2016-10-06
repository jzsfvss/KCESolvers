function ic = IC_sSweepCE(x, n)

global l
global u0

ic = [ u0(1)*MyHeaviside(x-l), u0(2)*Density(x/l, n), 0 ];