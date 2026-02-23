function save_tx_for_receiver(filename, tx_struct)
if ~(ischar(filename) || isstring(filename)), error('filename must be string'); end
if ~isfield(tx_struct,'tx_signal') || ~isfield(tx_struct,'X') || ~isfield(tx_struct,'params')
    error('tx_struct must contain tx_signal, X, and params');
end
tx_struct.tx_signal = tx_struct.tx_signal(:);
save(filename, '-struct', 'tx_struct');
fprintf('Saved %s (tx len=%d)\\n', filename, length(tx_struct.tx_signal));
end
