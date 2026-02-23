function run_all_tests()
addpath('tests');

test_dd_to_physical_mapping();
test_channel_scaling();
test_rvm_peak_location();

disp('All tests passed.');
end
