function den = D04_GaussianHeaviside(x)

den = MyHeaviside(x)*normpdf(x, 1.1, 0.4);
den = den*2/(1 - erf(-1.1/(0.4*sqrt(2)))); % Normalized due to Heaviside cut-off.