function plt = InvFigPlot(mth1, mth2, LpC0, LpC, u, t, rI, relI, ptit)

global mLpC0

% Linewidths and colors:
lnw1 = 4; % exact / L, T, C
lnw2 = 2; % numerical / L + C
lnw3 = 1; % error
%{
lnw1 = 5; % exact / L, T, C
lnw2 = 4; % numerical / L + C
lnw3 = 3; % error
%}
hue1 = 0.5; % exp. signal
hue2 = 1; % gen. signal
hue3 = 0.5; % error
colw = [1 1 1]; % white
colr = [1 0 0]; % red
coll = [0 0 1]; % blue
colg = [0 1 0]; % green
colb = [0 0 0]; % black
if (mth1 == 2); % MASKE
	col1 = colr;
	col2 = colg;
else % Standard or SimNECEEM
	col1 = colr;
	col2 = coll;
	col3 = colg;
end
if (mth1 ~= 2)
	Cind = 3;
	Cinds = [ 1, 3 ];
else
	Cind = 2;
	Cinds = [ 1, 2 ];
end

% Indices:
I = length(LpC0) - 1;
I1 = round(0.9*min(rI));
I2 = round(1.1*max(rI));
pI = I1:I2;
maxudet = max(LpC0);
sz = size(u);
J = sz(3)-1;
indJdet = J+1;

clf;
hold on;
axis([ t(I1), t(I2), 0, 1.1*maxudet ]);

% Experimental signal:
plot(t(pI), LpC0(pI), '-', 'Color', hue1*colb + (1-hue1)*colw, 'LineWidth', lnw1);

% Generated signal:
plot(t(pI), u(1, pI, indJdet), '-', 'Color', hue2*col1 + (1-hue2)*colw, 'LineWidth', lnw2);
plot(t(pI), u(2, pI, indJdet), '-', 'Color', hue2*col2 + (1-hue2)*colw, 'LineWidth', lnw2);
if (mth1 ~= 2)
	plot(t(pI), u(3, pI, indJdet), '-', 'Color', hue2*col3 + (1-hue2)*colw, 'LineWidth', lnw2);
end
plot(t(pI), u(1, pI, indJdet) + u(Cind, pI, indJdet), '-', 'Color', hue2*colb + (1-hue2)*colw, 'LineWidth', lnw3);

% Error:
rerrLC = 1.1*maxudet*abs(LpC0 - LpC)/mLpC0;
plot(t(pI), rerrLC(pI), ':', 'Color', hue3*colb + (1-hue3)*colw, 'LineWidth', lnw3);

% Delimiters:
l1 = vline(t(min(relI)), 'k:');
l2 = vline(t(max(relI)), 'k:');

% Legends:
if (mth1 == 2)
	%legend('L + C experimental', 'L generated', 'C generated', 'L + C generated', 'L + C rel. error', 't_min of inversion', 't_max of inversion');
	legend('L + C experimental', 'L fitted', 'C fitted', 'L + C fitted', 'L + C rel. error');
else
	%legend('L + C experimental', 'L generated', 'T generated', 'C generated', 'L + C generated', 'L + C rel. error', 't_min of inversion', 't_max of inversion');
	legend('L + C experimental', 'L fitted', 'T fitted', 'C fitted', 'L + C fitted', 'L + C rel. error');
end

% Other plot settings:
xlabel('Migration time (sec)');
ylabel('Concentration (mol / m^3)');
th = title(ptit);
P = get(th, 'Position');
set(th, 'Position', P.*[ 1, 1.015, 1 ]);

hold off;
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);

% Setting background to transparent:
set(gca, 'Color', 'None');
set(gca, 'TickDir', 'out');
glg = findobj(gcf,'Type','axes','Tag','legend');
set(glg, 'Color', 'None');
glg2 = findall(gcf, 'tag', 'legend');
set(glg2, 'location', 'northeastoutside', 'Color', 'None');

plt = 1;