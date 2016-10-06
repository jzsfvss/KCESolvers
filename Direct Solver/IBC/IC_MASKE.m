function ic = IC_MASKE(x, n)

global l
global u0

ic = u0*Density(x/l, n);