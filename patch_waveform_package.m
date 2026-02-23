function patch_waveform_package()
 % patch_waveform_package: copy top-level variables into params for receiver
 fn = 'waveform_package.mat';
 if ~isfile(fn), error('waveform_package.mat not found in current folder'); end
 S = load(fn);

 % ensure params exists
 if ~isfield(S,'params') || isempty(S.params)
     S.params = struct();
 end

 % move top-level targets -> params.targets
 if ~isfield(S.params,'targets') && isfield(S,'targets')
     S.params.targets = S.targets;
     fprintf('Copied top-level "targets" into params.targets\n');
 end

 % ensure grid sizes exist (infer from X)
 if (~isfield(S.params,'M') || ~isfield(S.params,'N')) && isfield(S,'X')
     [m,n] = size(S.X);
     S.params.M = m; S.params.N = n;
     fprintf('Inferred params.M=%d, params.N=%d from X\n',m,n);
 end
%% 

 % ensure Mcp present (default 0)
 if ~isfield(S.params,'Mcp'), S.params.Mcp = 0; end

 % Save back (use v7.3 for large matrices safety)
 save(fn,'-struct','S','-v7.3');
 fprintf('Saved patched %s\n', fn);
end
