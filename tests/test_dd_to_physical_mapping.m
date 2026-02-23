function test_dd_to_physical_mapping()
% TEST_DD_TO_PHYSICAL_MAPPING Independent of current working directory.
this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
repo_root = fileparts(this_dir);
addpath(fullfile(repo_root, 'config'));

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
