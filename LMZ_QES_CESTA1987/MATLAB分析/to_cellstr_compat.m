function c = to_cellstr_compat(x)
if iscell(x)
    c = x;
elseif ischar(x)
    c = cellstr(x);
elseif exist('isstring','builtin') && isstring(x)
    c = cellstr(x);
elseif isnumeric(x)
    c = arrayfun(@num2str,x,'UniformOutput',false);
else
    c = cellstr(x);
end
end
