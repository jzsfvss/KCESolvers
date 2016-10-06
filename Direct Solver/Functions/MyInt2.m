function y = MyInt2()

global mu
global eta

if (eta + mu == 0)
	y = 0;
else
	y = integral(@MyIntFun2, eta, -mu);
end