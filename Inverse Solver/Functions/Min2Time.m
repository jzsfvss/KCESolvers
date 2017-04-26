function tmend = Min2Time(mins)

rtend = clock;
rtend = rtend(4:5);
rttot = rtend(1)*60 + rtend(2) + round(mins);
hrs = floor(rttot/60);
mns = rttot - 60*hrs;
hrs = mod(hrs, 24);

if (hrs < 10)
	thrs = [ num2str(0), num2str(hrs) ];
else
	thrs = num2str(hrs);
end
if (mns < 10)
	tmns = [ num2str(0), num2str(mns) ];
else
	tmns = num2str(mns);
end

tmend = [ thrs, ':', tmns ];