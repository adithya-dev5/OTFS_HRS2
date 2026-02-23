function test_rvm_peak_location()
addpath('config'); addpath('core/sensing');
cfg = default_config();

RVM = zeros(64,64);
RVM(20,30) = 1;
peaks = find_rvm_peaks(RVM, 1, cfg);
assert(peaks(1,1)==20 && peaks(1,2)==30, 'RVM peak location mismatch.');
end
