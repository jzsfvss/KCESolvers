function [ rer, wrer, crer, lrer, lwrer ] = SolVerifierKCE(u, x, t, thr, threrr)

global k
global v
global D
global l
global u0

tind = length(t);
xind = length(x);
dt = abs(t(2)-t(1));
dx = abs(x(2)-x(1));
sz = (tind-2)*(xind-2);

M = zeros(tind-2, xind-2);
L = M;
T = M;
C = M;
L(:,:) = u(1, 2:(tind-1), 2:(xind-1));
T(:,:) = u(2, 2:(tind-1), 2:(xind-1));
C(:,:) = u(3, 2:(tind-1), 2:(xind-1));

dudt = (u(:, 3:tind, 2:(xind-1)) - u(:, 1:(tind-2), 2:(xind-1)))/(2*dt); % t = 2:(tind-1)
dudx = (u(:, 2:(tind-1), 3:xind) - u(:, 2:(tind-1), 1:(xind-2)))/(2*dx); % x = 2:(xind-1)
dudx2 = (u(:, 2:(tind-1), 3:xind) - 2*u(:, 2:(tind-1), 2:(xind-1)) + u(:, 2:(tind-1), 1:(xind-2)))/(dx^2); % x = 2:(xind-1)

lhsL = M;
lhsT = M;
lhsC = M;
lhsL(:,:) = dudt(1, :, :) + v(1)*dudx(1, :, :) - D(1)*dudx2(1, :, :);
lhsT(:,:) = dudt(2, :, :) + v(2)*dudx(2, :, :) - D(2)*dudx2(2, :, :);
lhsC(:,:) = dudt(3, :, :) + v(3)*dudx(3, :, :) - D(3)*dudx2(3, :, :);
rhsL = -k(1)*L.*T + k(2)*C;
rhsT = -k(1)*L.*T + k(2)*C;
rhsC = k(1)*L.*T - k(2)*C;

% Average absolute L2 error:
sz2 = size(lhsL);
J1 = 2;
J2 = sz2(2)-1;
sz = sz2(1)*(J2-J1+1);
aerL = sqrt(sum(sum((lhsL(:,J1:J2) - rhsL(:,J1:J2)).^2))/sz);
aerT = sqrt(sum(sum((lhsT(:,J1:J2) - rhsT(:,J1:J2)).^2))/sz);
aerC = sqrt(sum(sum((lhsC(:,J1:J2) - rhsC(:,J1:J2)).^2))/sz);
aer = [ aerL, aerT, aerC ];

% Average relative error:
mlhsL = max(max(abs(lhsL(:,J1:J2))));
mlhsT = max(max(abs(lhsT(:,J1:J2))));
mlhsC = max(max(abs(lhsC(:,J1:J2))));
mlhs = [ mlhsL, mlhsT, mlhsC ];
rer = 100*aer./mlhs;

% Average absolute weighted L2 error:
sz = size(L);
wL = zeros(sz(1), sz(2));
wT = wL;
wC = wL;
mL = max(max(L));
mT = max(max(T));
mC = max(max(C));
for i = 1:sz(1)
for j = 1:sz(2)
	wL(i,j) = (L(i,j) > thr*mL);
	wT(i,j) = (T(i,j) > thr*mT);
	wC(i,j) = (C(i,j) > thr*mC);
end
end
awerL = sqrt(sum(sum(wL(:,J1:J2).*(lhsL(:,J1:J2) - rhsL(:,J1:J2)).^2))/sum(sum(wL(:,J1:J2))));
awerT = sqrt(sum(sum(wT(:,J1:J2).*(lhsT(:,J1:J2) - rhsT(:,J1:J2)).^2))/sum(sum(wT(:,J1:J2))));
awerC = sqrt(sum(sum(wC(:,J1:J2).*(lhsC(:,J1:J2) - rhsC(:,J1:J2)).^2))/sum(sum(wC(:,J1:J2))));
awer = [ awerL, awerT, awerC ];

% Average relative error:
wrer = 100*awer./mlhs;

% Local relative error:
sz = size(lhsL);
lrer = zeros(3, sz(1), sz(2));
lwrer = zeros(3, sz(1), sz(2));
for (i = 1:sz(1))
for (j = 1:sz(2))
if (lhsL(i,j) ~= 0)
	lrer(1, i, j) = abs(lhsL(i,j)-rhsL(i,j));
	lwrer(1, i, j) = lrer(1, i, j)*wL(i,j);
else
	lrer(1, i, j) = 0;
	lwrer(1, i, j) = 0;
end
if (lhsT(i,j) ~= 0)
	lrer(2, i, j) = abs(lhsT(i,j)-rhsT(i,j));
	lwrer(2, i, j) = lrer(2, i, j)*wT(i,j);
else
	lrer(2, i, j) = 0;
	lwrer(2, i, j) = 0;
end
if (lhsC(i,j) ~= 0)
	lrer(3, i, j) = abs(lhsC(i,j)-rhsC(i,j));
	lwrer(3, i, j) = lrer(3, i, j)*wC(i,j);
else
	lrer(3, i, j) = 0;
	lwrer(3, i, j) = 0;
end
end
end
for l = 1:3
	lrer(l, :, :) = 100*lrer(l, :, :)/mlhs(l);
	lwrer(l, :, :) = 100*lwrer(l, :, :)/mlhs(l);
end

% Counted relative errors:
crer1 = 100*length(find(lwrer(1, :, :) > threrr))/length(find(lwrer(1, :, :) > thr));
crer2 = 100*length(find(lwrer(2, :, :) > threrr))/length(find(lwrer(2, :, :) > thr));
crer3 = 100*length(find(lwrer(3, :, :) > threrr))/length(find(lwrer(3, :, :) > thr));
crer = [ crer1, crer2, crer3 ];