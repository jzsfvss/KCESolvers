function [ rer, wrer, crer, lrer, lwrer ] = SolVerifierMASKE(u, x, t, thr, threrr)

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
A = M;
C = M;
A(:,:) = u(1, 2:(tind-1), 2:(xind-1));
C(:,:) = u(2, 2:(tind-1), 2:(xind-1));

dudt = (u(:, 3:tind, 2:(xind-1)) - u(:, 1:(tind-2), 2:(xind-1)))/(2*dt); % t = 2:(tind-1)
dudx = (u(:, 2:(tind-1), 3:xind) - u(:, 2:(tind-1), 1:(xind-2)))/(2*dx); % x = 2:(xind-1)
dudx2 = (u(:, 2:(tind-1), 3:xind) - 2*u(:, 2:(tind-1), 2:(xind-1)) + u(:, 2:(tind-1), 1:(xind-2)))/(dx^2); % x = 2:(xind-1)

lhsA = M;
lhsC = M;
lhsA(:,:) = dudt(1, :, :) + v(1)*dudx(1, :, :) - D(1)*dudx2(1, :, :);
lhsC(:,:) = dudt(2, :, :) + v(2)*dudx(2, :, :) - D(2)*dudx2(2, :, :);
rhsA = -k(1)*A*u0(2) + k(2)*C;
rhsC = k(1)*A*u0(2) - k(2)*C;

% Average absolute L2 error:
sz2 = size(lhsA);
J1 = 2;
J2 = sz2(2)-1;
sz = sz2(1)*(J2-J1+1);
aerA = sqrt(sum(sum((lhsA(:,J1:J2) - rhsA(:,J1:J2)).^2))/sz);
aerC = sqrt(sum(sum((lhsC(:,J1:J2) - rhsC(:,J1:J2)).^2))/sz);
aer = [ aerA, aerC ];

% Average relative error:
mlhsA = max(max(abs(lhsA(:,J1:J2))));
mlhsC = max(max(abs(lhsC(:,J1:J2))));
mlhs = [ mlhsA, mlhsC ];
rer = 100*aer./mlhs;

% Average absolute weighted L2 error:
sz = size(A);
wA = zeros(sz(1), sz(2));
wC = wA;
mA = max(max(A));
mC = max(max(C));
for i = 1:sz(1)
for j = 1:sz(2)
	wA(i,j) = (A(i,j) > thr*mA);
	wC(i,j) = (C(i,j) > thr*mC);
end
end
awerA = sqrt(sum(sum(wA(:,J1:J2).*(lhsA(:,J1:J2) - rhsA(:,J1:J2)).^2))/sum(sum(wA(:,J1:J2))));
awerC = sqrt(sum(sum(wC(:,J1:J2).*(lhsC(:,J1:J2) - rhsC(:,J1:J2)).^2))/sum(sum(wC(:,J1:J2))));
awer = [ awerA, awerC ];

% Average relative error:
wrer = 100*awer./mlhs;

% Local relative error:
sz = size(lhsA);
lrer = zeros(2, sz(1), sz(2));
lwrer = zeros(2, sz(1), sz(2));
for (i = 1:sz(1))
for (j = 1:sz(2))
if (lhsA(i,j) ~= 0)
	lrer(1, i, j) = abs(lhsA(i,j)-rhsA(i,j));
	lwrer(1, i, j) = lrer(1, i, j)*wA(i,j);
else
	lrer(1, i, j) = 0;
	lwrer(1, i, j) = 0;
end
if (lhsC(i,j) ~= 0)
	lrer(2, i, j) = abs(lhsC(i,j)-rhsC(i,j));
	lwrer(2, i, j) = lrer(2, i, j)*wC(i,j);
else
	lrer(2, i, j) = 0;
	lwrer(2, i, j) = 0;
end
end
end
for l = 1:2
	lrer(l, :, :) = 100*lrer(l, :, :)/mlhs(l);
	lwrer(l, :, :) = 100*lwrer(l, :, :)/mlhs(l);
end

% Counted relative errors:
crer1 = 100*length(find(lwrer(1, :, :) > threrr))/length(find(lwrer(1, :, :) > thr));
crer2 = 100*length(find(lwrer(2, :, :) > threrr))/length(find(lwrer(2, :, :) > thr));
crer = [ crer1, crer2 ];