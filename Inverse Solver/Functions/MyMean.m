function m = MyMean(x, n)

m = (sum(x.^n)/length(x))^(1/n);