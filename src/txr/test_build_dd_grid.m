%% test_build_dd_grid.m
clear; close all; clc;
% parameters identical to Midterm snippet (0-based lp,kp)
M = 32; N = 32;
lp = floor(M/2); kp = floor(N/2);
lTau = 2; kNu = 2;
Xp = 1+0j;

% create dummy DataSymbols based on mask counting
[dd_L, dd_K] = meshgrid(0:M-1, 0:N-1);
maskPilot = (dd_L == lp) & (dd_K == kp);
maskGuard = ( abs(dd_L - lp) <= lTau ) | ( abs(dd_K - kp) <= 2*kNu );
maskGuard = maskGuard & ~maskPilot;
maskData = ~(maskGuard | maskPilot);
numData = nnz(maskData);
% FIX: Use 2 (or any value != Xp) so we don't confuse Data with Pilot
DataSymbols = 2 * ones(numData,1);
% call build_dd_grid
X = build_dd_grid(M,N,lp,kp,lTau,kNu,DataSymbols,Xp);

% check pilot
pilot_val = X(kp+1, lp+1); % MATLAB 1-based read
if abs(pilot_val - Xp) ~= 0
    fprintf('FAIL: pilot value not equal to Xp.\n');
else
    fprintf('PASS: pilot value OK.\n');
end

% check guard/data mapping counts
countPilot = nnz(X == Xp);
countData  = nnz(X ~= 0 & X ~= Xp);
countTotal = N*M;
countGuard = countTotal - countPilot - countData;

% expected counts
expPilot = 1;
expGuard = nnz(maskGuard);
expData  = nnz(maskData);

if countPilot==expPilot && countGuard==expGuard && countData==expData
    fprintf('PASS: mask counts OK (pilot=%d, guard=%d, data=%d).\n',countPilot,countGuard,countData);
else
    fprintf('FAIL: mask count mismatch. got pilot=%d guard=%d data=%d | exp pilot=%d guard=%d data=%d\n', ...
        countPilot,countGuard,countData, expPilot, expGuard, expData);
end

% Reproduce DD_Grid plot from Midterm (visual check)
figure('Name','DD Grid Test'); hold on;
[idxK, idxL] = find(maskData); plot(idxL-1, idxK-1,'o','MarkerFaceColor',[0.55,0.35,0.75],'MarkerEdgeColor','k','MarkerSize',6);
[idxK, idxL] = find(maskGuard); plot(idxL-1, idxK-1,'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',6);
[idxK, idxL] = find(maskPilot); plot(idxL-1, idxK-1,'s','MarkerFaceColor',[1,0.55,0],'MarkerEdgeColor','k','MarkerSize',8);
axis equal; xlim([-0.5 M-0.5]); ylim([-0.5 N-0.5]); set(gca,'XTick',0:M-1,'YTick',0:N-1); grid on; box on;
xlabel('Delay'); ylabel('Doppler'); title('Pilot and Data structure (DD Grid)');
legend({'Data','Guard','Pilot'});
