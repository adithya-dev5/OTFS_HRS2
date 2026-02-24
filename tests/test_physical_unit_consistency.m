function test_physical_unit_consistency()
addpath('config'); addpath('core/utils');
cfg = default_config();
d = validate_physical_units(cfg);

c = 3e8;
range_per_delay_bin = c / (2*cfg.M*d.Delta_f);
vel_per_doppler_bin = c / (2*cfg.fc*cfg.N*d.T);

li = 7.0; ki = 3.0;
R = li * range_per_delay_bin;
V = ki * vel_per_doppler_bin;

li_back = (2*R*cfg.M*d.Delta_f)/c;
ki_back = (2*V*cfg.fc*cfg.N*d.T)/c;

assert(abs(li-li_back) < 1e-10, 'DD->range mapping inconsistency.');
assert(abs(ki-ki_back) < 1e-10, 'Doppler->velocity mapping inconsistency.');
end
