function y = MyRand(x, r)
% Generates a random value in the interval (x - r*x, x + r*x).

r2 = -1 + 2*rand(1);
r3 = 1 + r*r2;
y = r3*x;