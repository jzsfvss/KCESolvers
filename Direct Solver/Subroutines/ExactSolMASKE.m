function [ u, t, x, st, sx ] = ExactSolMASKE(IC0, ptp, I, J)

global mIC
global k
global v
global l
global u0
global xmax
global tmax
global ptp

global ome
global meps
global phi
global mpsi

global mu
global eta
global sig

% Mesh:
dt2 = tmax/I;
dx = xmax/J;
t = 0:dt2:tmax;
x = 0:dx:xmax;
I = length(t)-1;
J = length(x)-1;
st = I;
sx = J;

% Initialization:
u = zeros(2, I, J);
mIC = IC0;
vA = v(1);
vC = v(2);
B = u0(2);
kon = k(1);
koff = k(2);
Kd = koff/kon;
ome = 2*sqrt(kon*koff*u0(2))/(v(1)-v(2));
meps = sqrt(kon*B/koff);
phi = ome*(meps-(1/meps))/2;
mpsi = ome*(vC*meps-(vA/meps))/2;

% Calculation:
t1 = tic;
for i = 1:(I+1)
for j = 1:(J+1)

ti = t(i);
xj = x(j);
mu = v(1)*ti-xj;
eta = xj-v(2)*ti;
sig = ti*mpsi-xj*phi;

ic = mIC(xj-vA*ti, ptp);
icA = ic(1);
ic = mIC(xj-vC*ti, ptp);
icC = ic(3);

TA1 = icA*exp(-B*koff*ti/Kd);
TA2 = MyInt1();
TA3 = MyInt2();
u(1, i, j) = TA1 + TA2 + TA3;

TC1 = icC*exp(-koff*ti);
TC2 = MyInt3();
TC3 = MyInt4();
u(2, i, j) = TC1 + TC2 + TC3;

end

if (i == 50)
	t2 = toc(t1)/50;
	disp(['Estimated runtime = ', num2str(ceil(I*t2)), ' secs = ', num2str(MyRound(I*t2/60, 2)), ' mins.']);
end

end