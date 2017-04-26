function f = MetricFun(n)

switch (n)
case 1
	f = @Metric_L2;
case 2
	f = @Metric_MX;
otherwise
	f = @Metric_ExpL2;
end