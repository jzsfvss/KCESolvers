function a = MovAve(v, n, m)

I = (n-m):(n+m);

a = sum(v(I))/(2*m + 1);