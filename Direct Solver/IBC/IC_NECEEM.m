function ic = IC_NECEEM(x, n)

global l
global u0

ic = u0*Density(x/l, n);