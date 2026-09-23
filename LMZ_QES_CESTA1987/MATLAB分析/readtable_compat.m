function T = readtable_compat(fileName)
% readtable wrapper that works across older/newer MATLAB releases.
if ~exist(fileName,'file')
    error('File not found: %s', fileName);
end
try
    T = readtable(fileName,'VariableNamingRule','preserve');
catch
    try
        opts = detectImportOptions(fileName);
        T = readtable(fileName,opts);
    catch
        T = readtable(fileName);
    end
end
end
