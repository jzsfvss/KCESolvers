function se = SolExporter(eopt, plt, mth1, fnm, pnm, I, sx, u, u1, u2)

if (exist([ pwd, '\Exports']) == 0)
	system('mkdir Exports');
end

if (eopt == 0) % Export the entire solution.

if (plt == 1)
	sz = size(u);
else
	sz = size(u2);
end
tmpmtx = zeros(sz(2), sz(3));

if (mth1 == 2)
	nm = 'MASKE';
	snm = [ 'L', 'C' ];
	itn = 2;
else
	nm = 'SimNECEEM';
	snm = [ 'L', 'T', 'C' ];
	itn = 3;
end

if (plt == 1) % Exact or numerical solution export.

for (n = 1:itn)
	tmpmtx(:,:) = u(n, :, :);
	xlswrite([ '.\Exports\SolDS_', fnm, '_', snm(n) ], tmpmtx);
end
disp('Solution saved to files:');
disp([ pwd, '\Exports\SolDS_', fnm, '_*' ]);

else % Exact and numerical solution export.

fnm1 = [ 'SolDS_', nm, '_Exact_', pnm, '_', num2str(I), 'x', num2str(1) ];
for (n = 1:itn)
	tmpmtx(:, 1) = u1(n, :, 2);
	xlswrite([ '.\Exports\', fnm1, '_', snm(n) ], tmpmtx(:, 1));
end
disp('Exact solution saved to files:');
disp([ pwd, '\Exports\', fnm1, '_*' ]);

fnm2 = [ 'SolDS_', nm, '_Numerical_', pnm, '_', num2str(I), 'x', num2str(sx) ];
for (n = 1:itn)
	tmpmtx(:,:) = u2(n, :, :);
	xlswrite([ '.\Exports\', fnm2, '_', snm(n) ], tmpmtx);
end
disp('Numerical solution saved to files:');
disp([ pwd, '\Exports\', fnm2, '_*' ]);

end

else % Export the solution at the detector only.

if (mth1 == 2)
	nm = 'MASKE';
	snm = [ 'L', 'C' ];
	itn = 2;
else
	nm = 'SimNECEEM';
	snm = [ 'L', 'T', 'C' ];
	itn = 3;
end
snm2 = 'L+C';

if (plt == 1) % Exact or numerical solution export.

sz = size(u);
tmpmtx = zeros(sz(2), itn);
tmpmtx(:,:) = permute(u(1:itn, :, sz(3)), [ 2, 1, 3 ]);
xlswrite([ '.\Exports\SolDS_', fnm, '_', snm ], tmpmtx);
tmpmtx = zeros(sz(2), 1);
tmpmtx(:,:) = permute(u(1, :, sz(3)) + u(itn, :, sz(3)), [ 2, 1, 3 ]);
xlswrite([ '.\Exports\SolDS_', fnm, '_', snm2 ], tmpmtx);
disp('Solution saved to files:');
disp([ pwd, '\Exports\SolDS_', fnm, '_', snm, '.*' ]);
disp([ pwd, '\Exports\SolDS_', fnm, '_', snm2, '.*' ]);

else % Exact and numerical solution export.

fnm1 = [ 'SolDS_', nm, '_Exact_', pnm, '_', num2str(I), 'x', num2str(1) ];
sz = size(u1);
tmpmtx = zeros(sz(2), itn);
tmpmtx(:,:) = permute(u1(1:itn, :, sz(3)), [ 2, 1, 3 ]);
xlswrite([ '.\Exports\', fnm1, '_', snm ], tmpmtx);
tmpmtx = zeros(sz(2), 1);
tmpmtx(:,:) = permute(u1(1, :, sz(3)) + u1(itn, :, sz(3)), [ 2, 1, 3 ]);
xlswrite([ '.\Exports\', fnm1, '_', snm2 ], tmpmtx);
disp('Exact solution saved to files:');
disp([ pwd, '\Exports\', fnm1, '_', snm, '.*' ]);
disp([ pwd, '\Exports\', fnm1, '_', snm2, '.*' ]);

fnm2 = [ 'SolDS_', nm, '_Numerical_', pnm, '_', num2str(I), 'x', num2str(sx) ];
sz = size(u2);
tmpmtx = zeros(sz(2), itn);
tmpmtx(:,:) = permute(u2(1:itn, :, sz(3)), [ 2, 1, 3 ]);
xlswrite([ '.\Exports\', fnm2, '_', snm ], tmpmtx);
tmpmtx = zeros(sz(2), 1);
tmpmtx(:,:) = permute(u2(1, :, sz(3)) + u2(itn, :, sz(3)), [ 2, 1, 3 ]);
xlswrite([ '.\Exports\', fnm2, '_', snm2 ], tmpmtx);
disp('Numerical solution saved to files:');
disp([ pwd, '\Exports\', fnm2, '_', snm, '.*' ]);
disp([ pwd, '\Exports\', fnm2, '_', snm2, '.*' ]);

end

end