function ic = IC_ECEEM(x, n)

global l
global u0

ic = [ u0(1)*MyHeaviside(x-l), u0(2)*Density(x/l, n), u0(3)*MyHeaviside(x-l) ];