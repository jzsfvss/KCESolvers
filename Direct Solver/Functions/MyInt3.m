function y = MyInt3(ti, xj)

global mu
global eta

if (eta + mu == 0)
	y = 0;
else
	y = integral(@MyIntFun3, eta, -mu);
end