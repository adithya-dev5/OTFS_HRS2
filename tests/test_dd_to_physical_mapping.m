function test_dd_to_physical_mapping()
addpath('config');
cfg = default_config();

li = 5.3;
ki = 4.2;
tau = li / cfg.fs;
nu = ki * (cfg.fs/(cfg.M*cfg.N));

li_back = tau * cfg.fs;
ki_back = nu / (cfg.fs/(cfg.M*cfg.N));

assert(abs(li - li_back) < 1e-10, 'Delay index mapping mismatch.');
assert(abs(ki - ki_back) < 1e-10, 'Doppler index mapping mismatch.');
end
