function psd = read_oakridge_psd(fileName)
% Read the Oak Ridge UO2F2 initial PSD without requiring newer MATLAB APIs.
% Missing input files are reported as a warning and return an empty PSD.

cfg = CESTA1987_config();
psd = empty_psd();

if nargin < 1 || isempty(fileName)
    candidates = { ...
        fullfile(cfg.projectDir,'input','oakridge_psd.csv'), ...
        fullfile(cfg.projectDir,'oakridge_psd.csv'), ...
        fullfile(cfg.matlabDir,'input','oakridge_psd.csv'), ...
        fullfile(cfg.matlabDir,'oakridge_psd.csv'), ...
        fullfile(pwd,'input','oakridge_psd.csv'), ...
        fullfile(pwd,'oakridge_psd.csv')};
    fileName = '';
    for k = 1:numel(candidates)
        if exist(candidates{k},'file')
            fileName = candidates{k};
            break;
        end
    end
end

if isempty(fileName) || ~exist(fileName,'file')
    warning('Oak Ridge PSD not found. Checked input/oakridge_psd.csv and oakridge_psd.csv.');
    return;
end

try
    [diameter_um,values,meta] = read_psd_csv(fileName);
    [diameter_um,mass_fraction,number_fraction] = convert_distribution( ...
        diameter_um,values,meta);

    psd.exp.diameter_um = diameter_um(:);
    psd.exp.mass_fraction = mass_fraction(:);
    psd.exp.number_fraction = number_fraction(:);
    psd.exp.d50_um = mass_median_diameter(diameter_um,mass_fraction);
    psd.exp.sigma_g = mass_geometric_sigma(diameter_um,mass_fraction);
    psd.exp.diameter_type = meta.diameter_type;
    psd.exp.distribution_type = meta.distribution_type;
    psd.exp.measurement_time_s = meta.measurement_time_s;
    psd.exp.particle_density_kg_m3 = meta.particle_density_kg_m3;
    psd.exp.reference_density_kg_m3 = meta.reference_density_kg_m3;
    psd.exp.shape_factor = meta.shape_factor;
    psd.exp.source_file = fileName;
catch ME
    warning('Could not read Oak Ridge PSD %s: %s',fileName,ME.message);
    psd = empty_psd();
end
end

function psd = empty_psd()
psd = struct();
psd.exp = struct( ...
    'diameter_um',[], ...
    'mass_fraction',[], ...
    'number_fraction',[], ...
    'd50_um',NaN, ...
    'sigma_g',NaN, ...
    'diameter_type','undefined', ...
    'distribution_type','mass_fraction', ...
    'measurement_time_s',NaN, ...
    'particle_density_kg_m3',5000.0, ...
    'reference_density_kg_m3',1000.0, ...
    'shape_factor',1.0, ...
    'source_file','');
end

function [diameter_um,values,meta] = read_psd_csv(fileName)
fid = fopen(fileName,'r');
if fid < 0
    error('Unable to open file.');
end
cleanup = onCleanup(@() fclose(fid));

diameter_um = [];
values = [];
meta = struct( ...
    'diameter_type','undefined', ...
    'distribution_type','mass_fraction', ...
    'measurement_time_s',NaN, ...
    'particle_density_kg_m3',5000.0, ...
    'reference_density_kg_m3',1000.0, ...
    'shape_factor',1.0);

line = fgetl(fid);
while ischar(line)
    text = strtrim(line);
    if ~isempty(text)
        if text(1) == '#'
            meta = parse_directive(text,meta);
        else
            parts = strsplit(text,',');
            if numel(parts) >= 2
                dText = strtrim(parts{1});
                vText = strtrim(parts{2});
                if strcmpi(dText,'diameter_um') || strcmpi(dText,'diameter')
                    line = fgetl(fid);
                    continue;
                end
                d = str2double(dText);
                v = str2double(vText);
                if isfinite(d) && isfinite(v)
                    diameter_um(end+1,1) = d; %#ok<AGROW>
                    values(end+1,1) = v; %#ok<AGROW>
                end
            end
        end
    end
    line = fgetl(fid);
end

if isempty(diameter_um) || numel(diameter_um) ~= numel(values)
    error('No valid diameter_um,<distribution> rows found.');
end
if any(~isfinite(diameter_um)) || any(diameter_um <= 0)
    error('All PSD diameters must be finite and > 0.');
end
if any(~isfinite(values)) || any(values < 0)
    error('All PSD distribution values must be finite and >= 0.');
end
if ~any(values > 0)
    error('PSD distribution contains no positive values.');
end
end

function meta = parse_directive(text,meta)
eq = strfind(text,'=');
if isempty(eq)
    return;
end
key = lower(strtrim(text(2:eq(1)-1)));
value = strtrim(text(eq(1)+1:end));
switch key
    case 'diameter_type'
        meta.diameter_type = lower(value);
    case 'distribution_type'
        meta.distribution_type = lower(value);
    case 'measurement_time_s'
        meta.measurement_time_s = str2double(value);
    case 'particle_density'
        meta.particle_density_kg_m3 = str2double(value);
    case 'reference_density'
        meta.reference_density_kg_m3 = str2double(value);
    case 'shape_factor'
        meta.shape_factor = str2double(value);
end
end

function [diameter_um,mass_fraction,number_fraction] = convert_distribution( ...
    diameter_um,values,meta)
diameter_um = diameter_um(:);
values = values(:);

if strcmpi(meta.diameter_type,'aerodynamic')
    if ~isfinite(meta.particle_density_kg_m3) || meta.particle_density_kg_m3 <= 0 ...
            || ~isfinite(meta.reference_density_kg_m3) || meta.reference_density_kg_m3 <= 0 ...
            || ~isfinite(meta.shape_factor) || meta.shape_factor <= 0
        error('Invalid aerodynamic-to-geometric conversion parameters.');
    end
    diameter_um = diameter_um .* sqrt(meta.reference_density_kg_m3 ...
        ./ meta.particle_density_kg_m3) ./ sqrt(meta.shape_factor);
elseif ~strcmpi(meta.diameter_type,'geometric') ...
        && ~strcmpi(meta.diameter_type,'undefined')
    error('Unsupported diameter_type: %s',meta.diameter_type);
end

if any(diff(diameter_um) <= 0)
    error('PSD diameters must be strictly increasing.');
end

switch lower(meta.distribution_type)
    case 'mass_fraction'
        mass = values;
    case 'mass_density'
        width = log_interval_width(diameter_um);
        mass = values .* width;
    case 'number_density'
        width = log_interval_width(diameter_um);
        number = values .* width;
        diameter_m = diameter_um * 1.0e-6;
        particle_mass = meta.particle_density_kg_m3 * pi / 6.0 ...
            .* diameter_m.^3;
        mass = number .* particle_mass;
    otherwise
        error('Unsupported distribution_type: %s',meta.distribution_type);
end

massSum = sum(mass);
if ~isfinite(massSum) || massSum <= 0
    error('Converted PSD has zero mass.');
end
mass_fraction = mass / massSum;

diameter_m = diameter_um * 1.0e-6;
particle_mass = meta.particle_density_kg_m3 * pi / 6.0 .* diameter_m.^3;
number_raw = mass_fraction ./ particle_mass;
numberSum = sum(number_raw);
if isfinite(numberSum) && numberSum > 0
    number_fraction = number_raw / numberSum;
else
    number_fraction = zeros(size(mass_fraction));
end
end

function width = log_interval_width(diameter_um)
logDiameter = log(diameter_um(:));
n = numel(logDiameter);
if n < 2
    error('Density PSD integration requires at least two diameter points.');
end
width = zeros(n,1);
width(1) = logDiameter(2) - logDiameter(1);
for k = 2:n-1
    width(k) = 0.5 * (logDiameter(k+1) - logDiameter(k-1));
end
width(n) = logDiameter(n) - logDiameter(n-1);
if any(~isfinite(width)) || any(width <= 0)
    error('Invalid log-diameter integration width.');
end
end

function d50 = mass_median_diameter(diameter_um,mass_fraction)
d50 = NaN;
if isempty(diameter_um) || numel(diameter_um) ~= numel(mass_fraction)
    return;
end
p = mass_fraction(:) / sum(mass_fraction);
c = cumsum(p);
idx = find(c >= 0.5,1,'first');
if isempty(idx)
    return;
end
if idx == 1
    d50 = diameter_um(1);
    return;
end
c0 = c(idx-1);
if p(idx) <= 0
    d50 = diameter_um(idx);
else
    frac = (0.5 - c0) / p(idx);
    d50 = diameter_um(idx-1) + frac * (diameter_um(idx)-diameter_um(idx-1));
end
end

function sigma = mass_geometric_sigma(diameter_um,mass_fraction)
sigma = NaN;
q = isfinite(diameter_um) & isfinite(mass_fraction) ...
    & diameter_um > 0 & mass_fraction > 0;
if ~any(q)
    return;
end
p = mass_fraction(q) / sum(mass_fraction(q));
meanLog = sum(p .* log(diameter_um(q)));
sigma = exp(sqrt(sum(p .* (log(diameter_um(q))-meanLog).^2)));
end
