function [ lm, lmi, lcl ] = MyFindPeaks(v, rtol)

if (isempty(v))

lm = [];
lmi = [];
lcl = 0;

else

%[ lm, lmi ] = findpeaks(v, rtol);
if (length(v) > 2)
	[ lm, lmi ] = findpeaks(v);
	lcl = 1;
else
	lm = [];
	lcl = 0;
end

if (isempty(lm)) % No local peaks, so find the global one.
	[ lm, lmi ] = max(v);
	lcl = 0;
end % if 2

end % if 1