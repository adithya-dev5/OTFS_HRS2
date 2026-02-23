function Heff = build_heff(cfg)
% BUILD_HEFF Canonical Heff builder using unified channel_model interface.
MN = cfg.M * cfg.N;
Heff = zeros(MN, MN);
col_idx = 0;

for m = 1:cfg.M
    for n = 1:cfg.N
        col_idx = col_idx + 1;
        X_imp = zeros(cfg.N, cfg.M);
        X_imp(n,m) = 1;
        tx_imp = otfs_modulate(X_imp, cfg);
        if cfg.Mcp > 0
            tx_imp = [tx_imp(end-cfg.Mcp+1:end); tx_imp];
        end
        rx_imp = channel_model(tx_imp, cfg, cfg.channel.mode);
        Y_imp = otfs_demodulate(rx_imp, cfg);
        Heff(:, col_idx) = Y_imp(:);
    end
end
end
