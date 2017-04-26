function y = MyQuarticPlug(x0, c)

I1 = find((0 <= x0) & (x0 <= c));
I2 = find((c < x0) && (x0 <= 2*c));

x = 0*x0;
x(I1) = x0(I1);
x(I2) = 2*c - x0(I2);

y = (15/(4*(c^5)))*(x.^3).*((4*c/3)-x);