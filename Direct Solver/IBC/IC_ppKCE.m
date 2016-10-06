function ic = IC_ppKCE(x, n)

global l
global l2
global u0

ic = [ u0(1)*Density((x-l)/l2, n), u0(2)*Density(x/l2, n), 0 ];