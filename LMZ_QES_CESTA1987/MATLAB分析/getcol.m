function v = getcol(T, names, required)
% Case-insensitive table-column lookup with several fallback names.
if nargin < 3, required = true; end
if ischar(names), names = {names}; end
vn = T.Properties.VariableNames;
idx = [];
for k = 1:numel(names)
    j = find(strcmpi(vn,names{k}),1);
    if ~isempty(j), idx = j; break; end
end
if isempty(idx)
    if required
        error('Missing required column. Tried: %s', strjoin(names,', '));
    else
        v = [];
        return;
    end
end
v = T.(vn{idx});
end
