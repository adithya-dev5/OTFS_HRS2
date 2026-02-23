function Xhat = mmse_equalize(Y, Heff, noise_var)
% MMSE equalization (dense)
NM = numel(Y);
yvec = Y(:);

% regularization
A = Heff' * Heff + noise_var * eye(NM);
b = Heff' * yvec;

xhat = A \ b;
Xhat = reshape(xhat, size(Y));
end
