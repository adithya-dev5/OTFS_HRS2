%% OTFS Pilot + Guard + Data structure (Delay–Doppler grid)

% ---------------- Parameters (edit these) ----------------
M    = 32;          % # delay bins  (l = 0..M-1)   -> x-axis
N    = 32;          % # Doppler bins (k = 0..N-1)  -> y-axis
lp   = floor(M/2);  % pilot delay index (0-based)
kp   = floor(N/2);  % pilot Doppler index (0-based)
lTau = 2;           % guard half-width in delay (|l-lp| <= lTau)
kNu  = 2;           % guard half-width in Doppler (|k-kp| <= 2*kNu)

% ---------------- Build masks (0-based logical math) ----------------
% Coordinates for each cell
[L, K] = meshgrid(0:M-1, 0:N-1);

maskPilot = (L == lp) & (K == kp);
maskGuard = ((abs(L - lp) <= lTau) | (abs(K - kp) <= 2*kNu)) & ~maskPilot; % UNION (OR)
maskData  = ~(maskGuard | maskPilot);

% Optional: numeric grid if you want it returned/saved
% 0 = Guard, 1 = Data, 2 = Pilot
grids = zeros(N, M);
grids(maskData)  = 1;
grids(maskPilot) = 2;

% ---------------- Plot (markers like the paper) ----------------
figure; hold on;

% Plot data (purple circles)
[idxK, idxL] = find(maskData);
plot(idxL-1, idxK-1, 'o', 'MarkerFaceColor', [0.55 0.35 0.75], ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 6);

% Plot guard (white circles)
[idxK, idxL] = find(maskGuard);
plot(idxL-1, idxK-1, 'o', 'MarkerFaceColor', 'w', ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 6);

% Plot pilot (orange square)
[idxK, idxL] = find(maskPilot);
plot(idxL-1, idxK-1, 's', 'MarkerFaceColor', [1.0 0.55 0.0], ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 8);

% ---------------- Axes formatting to match Fig. 2 ----------------
axis equal; xlim([-0.5, M-0.5]); ylim([-0.5, N-0.5]);
set(gca, 'XTick', 0:M-1, 'YTick', 0:N-1);
grid on; box on; set(gca,'GridAlpha',0.25,'LineWidth',1.0);
xlabel('Delay'); ylabel('Doppler');
title('Pilot and Data symbols structure at the transmitter');

% Emphasize lp±lTau and kp±2kNu ticks (like the labels in the paper)
xt = unique([0 M-1 lp-lTau lp lp+lTau]);
yt = unique([0 N-1 kp-2*kNu kp kp+2*kNu]);
set(gca,'XTick', xt, 'YTick', yt);

% Legend
plot(nan,nan,'s','MarkerFaceColor',[1.0 0.55 0.0],'MarkerEdgeColor','k','MarkerSize',8);
plot(nan,nan,'o','MarkerFaceColor',[0.55 0.35 0.75],'MarkerEdgeColor','k','MarkerSize',6);
plot(nan,nan,'o','MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',6);
legend({'Guard','Data','Pilot'}, 'Location','northeastoutside');

% ---------------- Sanity check (optional) ----------------
fprintf('Pilot @ (k=%d, l=%d)\n', kp, lp);
fprintf('Guard columns: l in [%d, %d]\n', max(lp-lTau,0), min(lp+lTau,M-1));
fprintf('Guard rows:    k in [%d, %d]\n', max(kp-2*kNu,0), min(kp+2*kNu,N-1));
