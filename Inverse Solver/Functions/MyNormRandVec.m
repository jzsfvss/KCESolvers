function v = MyNormRandVec(n, mu, sig)

v = zeros(1, n);

for i = 1:n
	v(i) = normrnd(mu, sig);
end