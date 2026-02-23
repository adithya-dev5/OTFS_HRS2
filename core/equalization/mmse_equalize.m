function Xhat = mmse_equalize(Y, Heff, noise_var, cfg)
% MMSE_EQUALIZE Canonical MMSE equalizer with optional iterative stub.
if nargin < 4, cfg = struct(); end
NM = numel(Y);
yvec = Y(:);
A = Heff' * Heff + noise_var * eye(NM);
b = Heff' * yvec;

solver = "direct";
if isfield(cfg,'equalization') && isfield(cfg.equalization,'solver')
    solver = string(cfg.equalization.solver);
end

switch solver
    case "direct"
        xhat = A \ b;
    case "iterative"
        % Optional non-invasive path; falls back on direct if pcg unavailable.
        try
            [xhat, flag] = pcg(A, b, 1e-6, 200);
            if flag ~= 0
                warning('Iterative solver did not converge (flag=%d). Falling back to direct.', flag);
                xhat = A \ b;
            end
        catch
            warning('Iterative solver unavailable. Falling back to direct.');
            xhat = A \ b;
        end
    otherwise
        error('Unknown equalization solver: %s', solver);
end

Xhat = reshape(xhat, size(Y));
end
