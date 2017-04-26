function plt2 = FigPlot(plt, mth1, mth2, u, u1, u2, t, x, t1, x1, t2, x2)

% Linewidths and colors:
%{
lnw1 = 5; % exact / L, T, C
lnw2 = 2; % numerical / L + C
lnw3 = 2; % error
%}
lnw1 = 3; % exact / L, T, C
lnw2 = 1.5; % numerical / L + C
lnw3 = 1.5; % error
hue1 = 0.4; % exact
hue2 = 1; % numerical
hue3 = 0.6; % error
colw = [1 1 1]; % white
colr = [1 0 0]; % red
coll = [0 0 1]; % blue
colg = [0 1 0]; % green
colb = [0 0 0]; % black
if (mth1 == 2); % MASKE
	col1 = colr;
	col2 = colg;
else % KCE or SimNECEEM
	col1 = colr;
	col2 = coll;
	col3 = colg;
end

if (plt) % Single solution plot:

maxu = max(max(max(u)));
minu = min(min(min(u)));
sz = size(u);
indI = 1:sz(2);
indJ = 1:sz(3);

if ((mth1 == 2) && (mth2 == 2)) % MASKE Numerical.
	indJ = indJ(2:length(indJ));	
end
indJdet = indJ(length(indJ));
maxudet = max(max(max(u(:, :, indJdet))));
posind = find(u(1, :, indJdet) + u(2, :, indJdet) > 0.05*maxudet);
tmin0ind = max(1, min(posind) - ceil(20*sz(2)/1000));
tmax0ind = min(length(t), max(posind) + ceil(20*sz(2)/1000));
tmin0 = t(tmin0ind);
tmax0 = t(tmax0ind);
clf;
hold on;
axis([ tmin0, tmax0, 0, 1.1*maxudet ]);

if (mth1 == 2) % MASKE
	plot(t(indI), u(1, indI, indJdet), '-', 'Color', col1, 'LineWidth', lnw1);
	plot(t(indI), u(2, indI, indJdet), '-', 'Color', col2, 'LineWidth', lnw1);
	plot(t(indI), u(1, indI, indJdet) + u(2, indI, indJdet), '-', 'Color', colb, 'LineWidth', lnw2);
	legend('L', 'C', 'L + C');
else % KCE or SimNECEEM
	plot(t(indI), u(1, indI, indJdet), '-', 'Color', col1, 'LineWidth', lnw1);
	plot(t(indI), u(2, indI, indJdet), '-', 'Color', col2, 'LineWidth', lnw1);
	plot(t(indI), u(3, indI, indJdet), '-', 'Color', col3, 'LineWidth', lnw1);
	plot(t(indI), u(1, indI, indJdet) + u(3, indI, indJdet), '-', 'Color', colb, 'LineWidth', lnw2);
	legend('L', 'T', 'C', 'L + C');
end

xlabel('Migration time (sec)');
ylabel('Concentration (mol / m^3)');
hold off;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

else % Both solutions plot:

maxu = max(max(max(u2)));
minu = min(min(min(u2)));
sz = size(u2);
indI = 1:sz(2);
indJ = 1:sz(3);

maxudet1 = max(max(max(u1(:, :, 2))));
posind1 = find(u1(1, :, 2) + u1(2, :, 2) > 0.05*maxudet1);
tmin0ind1 = max(1, min(posind1) - ceil(20*sz(2)/1000));
tmax0ind1 = min(length(t), max(posind1) + ceil(20*sz(2)/1000));
tmin01 = t1(tmin0ind1);
tmax01 = t1(tmax0ind1);

indJdet = indJ(length(indJ));
maxudet2 = max(max(max(u2(:, :,indJdet))));
posind2 = find(u2(1, :, indJdet) + u2(2, :, indJdet) > 0.05*maxudet2);
tmin0ind2 = max(1, min(posind2) - ceil(20*sz(2)/1000));
tmax0ind2 = min(length(t), max(posind2) + ceil(20*sz(2)/1000));
tmin02 = t2(tmin0ind2);
tmax02 = t2(tmax0ind2);

maxudet = max(maxudet1, maxudet2);
tmin0 = min(tmin01, tmin02);
tmax0 = max(tmax01, tmax02);

clf;
hold on;
axis([ tmin0, tmax0, 0, 1.1*maxudet ]);

% Exact:
plot(t1(indI), u1(1, indI, 2), '-', 'Color', hue1*col1 + (1-hue1)*colw, 'LineWidth', lnw1);
plot(t1(indI), u1(2, indI, 2), '-', 'Color', hue1*col2 + (1-hue1)*colw, 'LineWidth', lnw1);
if (mth1 == 2)
	plot(t1(indI), u1(1, indI, 2) + u1(2, indI, 2), '-', 'Color', hue1*colb + (1-hue1)*colw, 'LineWidth', lnw1);
else
	plot(t1(indI), u1(3, indI, 2), '-', 'Color', hue1*col3 + (1-hue1)*colw, 'LineWidth', lnw1);
end

% Numerical:
plot(t2(indI), u2(1, indI, indJdet), '-', 'Color', hue2*col1 + (1-hue2)*colw, 'LineWidth', lnw2);
plot(t2(indI), u2(2, indI, indJdet), '-', 'Color', hue2*col2 + (1-hue2)*colw, 'LineWidth', lnw2);
if (mth1 == 2)
	plot(t2(indI), u2(1, indI, indJdet) + u2(2, indI, indJdet), '-', 'Color', hue2*colb + (1-hue2)*colw, 'LineWidth', lnw2);
else
	plot(t2(indI), u2(3, indI, indJdet), '-', 'Color', hue2*col3 + (1-hue2)*colw, 'LineWidth', lnw2);
end

% Error:
if (mth1 == 2)

u2intL = interp1(t2(indI), u2(1, indI, indJdet), t1(indI), 'pchip', 'extrap');
u2intC = interp1(t2(indI), u2(2, indI, indJdet), t1(indI), 'pchip', 'extrap');
u2intLC = interp1(t2(indI), u2(1, indI, indJdet) + u2(2, indI, indJdet), t1(indI), 'pchip', 'extrap');

mxL = abs(max(max(u1(1, indI, 2))));
rerrL = 1.1*maxudet*abs(u1(1, indI, 2) - u2intL)/mxL;
mxC = abs(max(max(u1(2, indI, 2))));
rerrC = 1.1*maxudet*abs(u1(2, indI, 2) - u2intC)/mxC;
mxLC = abs(max(max(u1(1, indI, 2) + u1(2, indI, 2))));
rerrLC = 1.1*maxudet*abs((u1(1, indI, 2) + u1(2, indI, 2)) - u2intLC)/mxLC;

else

u2intL = interp1(t2(indI), u2(1, indI, indJdet), t1(indI), 'pchip', 'extrap');
u2intT = interp1(t2(indI), u2(2, indI, indJdet), t1(indI), 'pchip', 'extrap');
u2intC = interp1(t2(indI), u2(3, indI, indJdet), t1(indI), 'pchip', 'extrap');
u2intLC = interp1(t2(indI), u2(1, indI, indJdet) + u2(3, indI, indJdet), t1(indI), 'pchip', 'extrap');

mxL = abs(max(max(u1(1, indI, 2))));
rerrL = 1.1*maxudet*abs(u1(1, indI, 2) - u2intL)/mxL;
mxT = abs(max(max(u1(2, indI, 2))));
rerrT = 1.1*maxudet*abs(u1(2, indI, 2) - u2intT)/mxT;
mxC = abs(max(max(u1(3, indI, 2))));
rerrC = 1.1*maxudet*abs(u1(3, indI, 2) - u2intC)/mxC;
mxLC = abs(max(max(u1(1, indI, 2) + u1(3, indI, 2))));
rerrLC = 1.1*maxudet*abs((u1(1, indI, 2) + u1(3, indI, 2)) - u2intLC)/mxLC;

end

plot(t1(indI), rerrL, ':', 'Color', hue3*col1 + (1-hue3)*colw, 'LineWidth', lnw3);
if (mth1 == 2)
	plot(t1(indI), rerrC, ':', 'Color', hue3*col2 + (1-hue3)*colw, 'LineWidth', lnw3);
	plot(t1(indI), rerrLC, ':', 'Color', hue3*colb + (1-hue3)*colw, 'LineWidth', lnw3);
else
	plot(t1(indI), rerrT, ':', 'Color', hue3*col2 + (1-hue3)*colw, 'LineWidth', lnw3);
	plot(t1(indI), rerrC, ':', 'Color', hue3*col3 + (1-hue3)*colw, 'LineWidth', lnw3);
end

if (mth1 == 2)
	legend('L exact', 'C exact', 'L + C exact', 'L numerical', 'C numerical', 'L + C numerical', 'L rel. error', 'C rel. error', 'L + C rel. error');
else
	legend('L exact', 'T exact', 'C exact', 'L numerical', 'T numerical', 'C numerical', 'L rel. error', 'T rel. error', 'C rel. error');
end
xlabel('Migration time (sec)');
ylabel('Concentration (mol / m^3)');
hold off;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

end

% Setting background to transparent:
set(gca, 'Color', 'None');
set(gca, 'TickDir', 'out');
glg = findobj(gcf, 'Type', 'axes', 'Tag', 'legend');
set(glg, 'Color', 'None');
glg2 = findall(gcf, 'tag', 'legend');
set(glg2, 'location', 'northeastoutside', 'Color', 'None');

plt2 = 1;