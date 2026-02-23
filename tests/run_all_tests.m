function run_all_tests()
% RUN_ALL_TESTS Robust test runner from any working directory.
this_file = mfilename('fullpath');
if isempty(this_file)
    this_dir = pwd;
else
    this_dir = fileparts(this_file);
end
repo_root = fileparts(this_dir);

addpath(this_dir);
addpath(fullfile(repo_root, 'config'));
addpath(fullfile(repo_root, 'core', 'utils'));
addpath(fullfile(repo_root, 'core', 'channel'));
addpath(fullfile(repo_root, 'core', 'sensing'));

test_dd_to_physical_mapping();
test_channel_scaling();
test_rvm_peak_location();
test_physical_unit_consistency();

disp('All tests passed.');
end
