function y = MyInt1()

global mu
global eta

if (eta + mu == 0)
	y = 0;
else
	y = integral(@MyIntFun1, eta, -mu);
end