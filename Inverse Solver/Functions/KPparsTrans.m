function kpptr = KPparsTrans(kpp, cT, sT)

global l0
global Lini
global Tini

% Initialization:
fac = 10; % Factor for order of magnitude bound.
koff = kpp(1);
sL1 = kpp(3);
sL2 = kpp(4);
hL = kpp(5);
sC1 = kpp(8);
sC2 = kpp(9);
hC = kpp(10);
if ((kpp(6) > cT*fac) || (kpp(6) < cT/fac)) % Prevents runaway T plug center.
	kpp(6) = cT;
end
sT1 = sT(1);
sT2 = sT(2);

% Eq. concentrations:
Leq = (hL/l0)*sqrt(2*pi)*(sL1 + sL2)/2;
Ceq = (hC/l0)*sqrt(2*pi)*(sC1 + sC2)/2;
Teq = 1000*Tini - Ceq;

% Calculation:
sTa = (sT1 + sT2)/2;
hT = l0*Teq/(sqrt(2*pi)*sTa);
kon = koff*Ceq/(Leq*Teq);

% Assignment:
kpptr = abs([ kon, kpp(1:6), sT, hT, kpp(7:10) ]);