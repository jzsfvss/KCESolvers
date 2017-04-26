function [ u, t, x, st, sx ] = ExactSolSimNECEEM(IC0, ptp, I, J, errtol)

global k
global v
global D
global l
global u0
global xmax
global tmax

% Mesh:
t = linspace(0, tmax, I+1);
x = linspace(0, xmax, J+1);
I = length(t)-1;
J = length(x)-1;
st = I;
sx = J;

% Initialization of variables:
u = zeros(3, I+1, J+1);
vL = v(1);
vT = v(2);
vC = v(3);
DL = D(1);
DT = D(2);
DC = D(3);
L0 = u0(1);
T0 = u0(2);
C0 = u0(3);
koff = k(2);

% Initial condition:
for j = 1:(J+1)
	ic = IC0(x(j), ptp);
	u(1, 1, j) = ic(1);
	u(2, 1, j) = ic(2);
	u(3, 1, j) = ic(3);
end

if (ptp == 1) % Heaviside.

t1 = tic;
for i = 2:(I+1)
	ti = t(i);
	for j = 1:(J+1)
		xj = x(j);
		% Calculating L and T:
		for n = 1:2
			unijeq = (u0(n)/2)*(erf((l-xj+ti*v(n))/sqrt(4*ti*D(n))) - erf((ti*v(n)-xj)/sqrt(4*ti*D(n))));
			h = @(r,z) exp(-r*koff).*exp(-0.25*((z+(vC-v(n))*r-xj+v(n)*ti).^2)./(DC*r+D(n)*(ti-r)))./sqrt(4*pi*(DC*r+D(n)*(ti-r)));
			unijdis = (koff*C0)*integral2(h, 0, ti, 0, l, 'Method', 'iterated', 'AbsTol', 0, 'RelTol', errtol);
			u(n, i, j) = unijeq + unijdis;
		end
		% Calculating C:
		u(3, i, j) = (C0/2)*exp(-ti*koff)*(erf((l-xj+ti*vC)/sqrt(4*ti*DC)) - erf((ti*vC-xj)/sqrt(4*ti*DC)));
	end
	if (i == 50)
		t2 = toc(t1)/50;
		disp(['Estimated runtime = ', num2str(ceil(I*t2)), ' secs = ', num2str(MyRound(I*t2/60, 2)), ' mins.']);
	end
end

else % Gaussian and Shifted-Gaussian.

if (ptp == 2)
	m0 = 1.1;
	s0 = 0.4;
else
	m0 = 1.5;
	s0 = 0.4;
end

t1 = tic;
for i = 2:(I+1)
	ti = t(i);
	for j = 1:(J+1)
		xj = x(j);
		% Calculating L and T:
		for n = 1:2
			unijeq = l*u0(n)*normpdf(xj, m0*l + v(n)*ti, sqrt(((s0*l)^2) + 2*D(n)*ti));
			h = @(r) exp(-koff*r).*normpdf(xj, m0*l + v(n)*ti + (vC-v(n))*r, sqrt(((s0*l)^2) + 2*D(n)*ti + 2*(DC-D(n))*r));
			unijdis = koff*l*C0*integral(h, 0, ti, 'RelTol', errtol, 'AbsTol', 0);
			%unijdis = koff*l*C0*integral(h, 0, ti, 'RelTol', 1E-6);
			u(n, i, j) = unijeq + unijdis;
		end
		% Calculating C:
		u(3, i, j) = l*C0*exp(-koff*ti)*normpdf(xj, m0*l + vC*ti, sqrt(((s0*l)^2) + 2*DC*ti));
	end
	if (i == 50)
		t2 = toc(t1)/50;
		disp(['Estimated runtime = ', num2str(ceil(I*t2)), ' secs = ', num2str(MyRound(I*t2/60, 2)), ' mins.']);
	end
end

end