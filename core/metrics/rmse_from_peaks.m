function out = rmse_from_peaks(peaks, true_targets, cfg)
% RMSE_FROM_PEAKS RMSE in DD index space with canonical axis mapping.
assert(cfg.dd_axis == "doppler_row_delay_col", 'Unsupported dd_axis convention.');
out = struct('RMSE_range_idx', NaN, 'RMSE_doppler_idx', NaN);
if isempty(peaks) || isempty(true_targets), return; end

Lmin = min(size(peaks,1), numel(true_targets));
est_ki = peaks(1:Lmin,1);
est_li = peaks(1:Lmin,2);

true_li = nan(Lmin,1);
true_ki = nan(Lmin,1);
for i=1:Lmin
    tg = true_targets(i);
    if isfield(tg,'li'), true_li(i) = tg.li;
    elseif isfield(tg,'tau'), true_li(i) = round(tg.tau * cfg.fs); end

    if isfield(tg,'ki'), true_ki(i) = tg.ki;
    elseif isfield(tg,'nu'), true_ki(i) = round(tg.nu / (cfg.fs/(cfg.M*cfg.N))); end
end

vli = ~isnan(true_li);
vki = ~isnan(true_ki);
if any(vli), out.RMSE_range_idx = sqrt(mean((est_li(vli)-true_li(vli)).^2)); end
if any(vki), out.RMSE_doppler_idx = sqrt(mean((est_ki(vki)-true_ki(vki)).^2)); end
end
