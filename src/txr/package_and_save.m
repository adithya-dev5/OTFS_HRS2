%% package_and_save.m
cfg = struct('M',32,'N',32,'lp',floor(32/2)-1,'kp',floor(32/2)-1,'lTau',2,'kNu',2,'Xp',1+0j,'fs',15.36e6,'Mcp',64);
[X, tx, params, ~] = waveform_and_channel(cfg);
targets(1).tau = 1.2e-6; targets(1).nu = 3; targets(1).gain = 0.8*exp(1j*0.2);
targets(2).tau = 3.8e-6; targets(2).nu = -1.5; targets(2).gain = 0.4*exp(-1j*0.5);
waveform_package.X = X; waveform_package.tx_signal = tx; waveform_package.params = params;
waveform_package.targets = targets; waveform_package.Lfft = params.Lfft;
waveform_package.fvec = params.fvec; waveform_package.tvec = params.tvec;
save('waveform_package.mat','-struct','waveform_package');
fprintf('Saved waveform_package.mat into current folder.\\n');
