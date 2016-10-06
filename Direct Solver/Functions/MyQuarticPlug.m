function y = MyQuarticPlug(x0, c)

if ((0 <= x0) && (x0 <= c))
	x = x0;
elseif ((c < x0) && (x0 <= 2*c))
	x = 2*c-x0;
else
	x = 0;
end

y = (15/(4*(c^5)))*(x.^3).*((4*c/3)-x);