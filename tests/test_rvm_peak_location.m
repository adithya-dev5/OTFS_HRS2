function test_rvm_peak_location()
% TEST_RVM_PEAK_LOCATION Independent of current working directory.
this_file = mfilename('fullpath');
this_dir = fileparts(this_file);
repo_root = fileparts(this_dir);
addpath(fullfile(repo_root, 'config'));
addpath(fullfile(repo_root, 'core', 'sensing'));

cfg = default_config();

RVM = zeros(64,64);
RVM(20,30) = 1;
peaks = find_rvm_peaks(RVM, 1, cfg);
assert(peaks(1,1)==20 && peaks(1,2)==30, 'RVM peak location mismatch.');
end
