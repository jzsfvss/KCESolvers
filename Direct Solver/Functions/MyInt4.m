function y = MyInt4()

global mu
global eta

if (eta + mu == 0)
	y = 0;
else
	y = integral(@MyIntFun4, eta, -mu);
end