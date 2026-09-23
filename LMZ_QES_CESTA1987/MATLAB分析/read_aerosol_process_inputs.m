function inputs = read_aerosol_process_inputs()
% Read optional aerosol process diagnostics produced by the QES run.
% Missing files are reported but do not stop the remaining PSD analysis.

cfg = CESTA1987_config();
inputs = struct();
[inputs.psd,inputs.psdFile] = read_optional( ...
    cfg,'aerosol_psd.csv','aerosol PSD particle output');
[inputs.mass,inputs.massFile] = read_optional( ...
    cfg,'aerosol_mass_diagnostic.csv','aerosol mass diagnostic');
[inputs.deposition,inputs.depositionFile] = read_optional( ...
    cfg,'aerosol_deposition.csv','aerosol deposition inventory');
inputs.hasPsd = ~isempty(inputs.psd);
inputs.hasMass = ~isempty(inputs.mass);
inputs.hasDeposition = ~isempty(inputs.deposition);
end

function [T,fileName] = read_optional(cfg,baseName,description)
T = [];
fileName = find_output_file(cfg,baseName);
if isempty(fileName)
    warning('%s not found: %s',description,baseName);
    return;
end
try
    T = readtable_compat(fileName);
catch ME
    warning('Could not read %s (%s): %s',description,fileName,ME.message);
    T = [];
end
end

function fileName = find_output_file(cfg,baseName)
candidates = { ...
    fullfile(cfg.outputDir,baseName), ...
    fullfile(cfg.projectDir,baseName), ...
    fullfile(cfg.matlabDir,baseName), ...
    fullfile(pwd,baseName)};
fileName = '';
for k = 1:numel(candidates)
    if exist(candidates{k},'file')
        fileName = candidates{k};
        return;
    end
end
end
