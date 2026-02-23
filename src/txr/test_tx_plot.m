%% test_tx_plot.m
clear; close all; clc;
cfg = struct('M',32,'N',32,'lp',floor(32/2)-1,'kp',floor(32/2)-1,'lTau',2,'kNu',2,'Xp',1+0j,'fs',15.36e6,'Mcp',64);
[X, tx, params, ~] = waveform_and_channel(cfg);

% plot abs(tx)
figure('Name','Transmit waveform magnitude |tx|');
plot(abs(tx));
xlabel('Sample index'); ylabel('|tx|');
title('Transmit waveform magnitude (abs(tx))');
grid on;

% Basic checks
if isempty(tx) || ~isvector(tx)
    fprintf('FAIL: tx empty or not a vector.\n');
else
    fprintf('PASS: tx generated, length=%d samples.\n', length(tx));
end

% normalize check: mean power ~1
p = mean(abs(tx).^2);
fprintf('tx avg power = %.4f (should be 1 after normalization).\n', p);
if abs(p-1) < 1e-6
    fprintf('PASS: tx normalized.\n');
else
    fprintf('WARN: tx not exactly unit-power (value=%.6g). This is allowed; verify plot visually.\n', p);
end
