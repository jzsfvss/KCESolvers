function v = MyIntFun3(Y)

global mIC
global ome
global phi
global mu
global eta
global sig
global ptp

lY = length(Y);
v = zeros(1,lY);

for i = 1:lY

y = Y(i);
ic = mIC(y, ptp);
v(i) = abs(ome/2)*sqrt((mu+y)/(eta-y))*besseli(1, abs(ome)*sqrt((mu+y)*(eta-y)))*exp(sig+y*phi)*ic(3);
v(i) = real(v(i));

end