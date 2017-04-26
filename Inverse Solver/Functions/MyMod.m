function y = MyMod(x, m)

y = mod(x, m);

if (y == 0)
	y = m;
end