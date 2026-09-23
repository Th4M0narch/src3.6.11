function psd = read_qes_psd_sampler(fileName)
% Read only aerosol_sampler.csv: nine aerodynamic observation classes per
% collector/time. Diameters in the file are metres and use dry-mass weighting;
% optional number-weighted columns remain distinct. Never read parcel history.

psd = empty_psd();
if nargin < 1 || isempty(fileName)
    cfg = CESTA1987_config();
    fileName = cfg.aerosolSamplerCsv;
end
if ~exist(fileName,'file')
    warning('QES:PSD:MissingFile','QES aerosol sampler CSV not found: %s',fileName);
    return;
end
T = readtable_compat(fileName);
if isempty(T)
    warning('QES:PSD:EmptyFile','QES aerosol sampler CSV contains no rows: %s',fileName);
    return;
end
% A component-history or legacy wide table is a schema error, not a zero PSD.
time = required_numeric(T,'time');
ids = required_text(T,'collector_id');
bins = required_numeric(T,'bin_id');
mass = required_numeric(T,'mass_fraction');
number = required_numeric(T,'number_fraction');
dry = required_numeric(T,'dry_d');
wet = required_numeric(T,'wet_d');
aero = required_numeric(T,'aerodynamic_d');
x = required_numeric(T,'x');
y = required_numeric(T,'y');
z = required_numeric(T,'z');
w0 = required_numeric(T,'window_start');
w1 = required_numeric(T,'window_end');
duration = required_numeric(T,'sampling_time_s');
complete = required_numeric(T,'sample_complete');
units = required_text(T,'diameter_unit');
basis = required_text(T,'diameter_basis');
massBasis = required_text(T,'mass_basis');
if any(~strcmp(units,'m')) || any(~strcmp(basis,'aerodynamic')) ...
        || any(~strcmp(massBasis,'dry_UO2F2'))
    error('QES:PSD:Schema','Expected m, aerodynamic, and dry_UO2F2 metadata.');
end
values = [time bins mass number dry wet aero x y z w0 w1 duration complete];
if any(~isfinite(values(:))) || any(time < 0) || any(w0 < 0) ...
        || any(w1 <= w0) || any(duration < 0) || any(duration > w1-w0+1e-8) ...
        || any(~ismember(complete,[0 1])) || any(bins ~= fix(bins)) ...
        || any(bins < 1 | bins > 9) || any(mass < 0 | number < 0) ...
        || any(dry < 0 | wet < 0 | aero < 0) || any(cellfun(@isempty,ids))
    error('QES:PSD:Schema','Invalid collector metadata, fractions, or geometry.');
end
if any(complete == 1 & (time < w1-1e-8 | abs(duration-(w1-w0)) > 1e-8))
    error('QES:PSD:Schema','A completed sample must cover its full collection window.');
end
populated = mass > 0 | number > 0;
if any(populated & (dry <= 0 | wet <= 0 | aero <= 0)) ...
        || any((mass > 0) ~= (number > 0))
    error('QES:PSD:Schema','Populated classes require positive mass, number and diameters.');
end
distance = optional_numeric(T,'distance_m');
azimuth = optional_numeric(T,'azimuth_deg');
if isempty(distance), distance = hypot(x-550.0,y-100.0); end
if isempty(azimuth), azimuth = mod(atan2d(x-550.0,y-100.0),360); end
if any(~isfinite(distance)) || any(distance < 0) || any(~isfinite(azimuth))
    error('QES:PSD:Schema','Invalid CESTA distance or azimuth.');
end
% findgroups groups numeric times without lossy rounding or periodic sampling.
group = findgroups(ids,time);
groups = unique(group);
n = numel(groups);
q = psd.qes;
q.mass_fraction = zeros(n,9);
q.number_fraction = zeros(n,9);
q.dry_d_m = zeros(n,9);
q.wet_d_m = zeros(n,9);
q.aerodynamic_d_m = zeros(n,9);
q.dry_d_number_weighted_m = NaN(n,9);
q.wet_d_number_weighted_m = NaN(n,9);
q.aerodynamic_d_number_weighted_m = NaN(n,9);
d50 = optional_numeric(T,'d50_um');
sigma = optional_numeric(T,'sigma_g');
numberDiameters = {optional_numeric(T,'dry_d_number_weighted'), ...
    optional_numeric(T,'wet_d_number_weighted'), ...
    optional_numeric(T,'aerodynamic_d_number_weighted')};
numberFields = {'dry_d_number_weighted_m','wet_d_number_weighted_m', ...
    'aerodynamic_d_number_weighted_m'};
for k = 1:n
    rows = find(group == groups(k));
    [sortedBins,order] = sort(bins(rows));
    rows = rows(order);
    if numel(rows) ~= 9 || ~isequal(sortedBins(:),(1:9)')
        error('QES:PSD:Schema','Each collector/time must contain bins 1..9 exactly once.');
    end
    first = rows(1);
    metadata = [x(rows) y(rows) z(rows) distance(rows) azimuth(rows) ...
        w0(rows) w1(rows) duration(rows) complete(rows)];
    if any(any(abs(metadata-repmat(metadata(1,:),9,1)) > 1e-10))
        error('QES:PSD:Schema','Metadata differ between bins of one collector/time.');
    end
    validate_fraction(mass(rows));
    validate_fraction(number(rows));
    q.time_s(k,1) = time(first);
    q.sampler_id{k,1} = ids{first};
    q.distance_m(k,1) = distance(first);
    q.height_m(k,1) = z(first);
    q.azimuth_deg(k,1) = mod(azimuth(first),360);
    q.window_start(k,1) = w0(first);
    q.window_end(k,1) = w1(first);
    q.sampling_time_s(k,1) = duration(first);
    q.sample_complete(k,1) = logical(complete(first));
    q.mass_fraction(k,:) = mass(rows)';
    q.number_fraction(k,:) = number(rows)';
    q.dry_d_m(k,:) = dry(rows)';
    q.wet_d_m(k,:) = wet(rows)';
    q.aerodynamic_d_m(k,:) = aero(rows)';
    for d = 1:3
        values = numberDiameters{d};
        if ~isempty(values)
            if any(~isfinite(values(rows))) || any(values(rows) < 0) ...
                    || any(number(rows) > 0 & values(rows) <= 0)
                error('QES:PSD:Schema','Invalid number-weighted diameter.');
            end
            q.(numberFields{d})(k,:) = values(rows)';
        end
    end
    if isempty(d50)
        index = find(cumsum(mass(rows)) >= 0.5,1);
        q.d50_um(k,1) = NaN;
        if ~isempty(index), q.d50_um(k,1) = aero(rows(index))*1e6; end
    else
        consistent_statistic(d50(rows));
        q.d50_um(k,1) = d50(first);
    end
    if isempty(sigma)
        active = mass(rows) > 0;
        weights = mass(rows(active));
        logs = log(aero(rows(active)));
        q.sigma_g(k,1) = NaN;
        if ~isempty(weights)
            meanLog = sum(weights.*logs)/sum(weights);
            q.sigma_g(k,1) = exp(sqrt(sum(weights.*(logs-meanLog).^2)/sum(weights)));
        end
    else
        consistent_statistic(sigma(rows));
        q.sigma_g(k,1) = sigma(first);
    end
end
q.source_file = fileName;
q.diameter_basis = 'aerodynamic';
q.mass_basis = 'dry_UO2F2';
psd.qes = q;
end

function psd = empty_psd()
psd.qes = struct('time_s',[],'distance_m',[],'height_m',[], ...
    'azimuth_deg',[],'window_start',[],'window_end',[], ...
    'sampling_time_s',[],'sample_complete',[], ...
    'mass_fraction',zeros(0,9),'number_fraction',zeros(0,9), ...
    'dry_d_m',zeros(0,9),'wet_d_m',zeros(0,9),'aerodynamic_d_m',zeros(0,9), ...
    'd50_um',[],'sigma_g',[],'sampler_id',{{}}, ...
    'bin_diameter_um',[0.2 sqrt(.4*.7) sqrt(.7*1.1) sqrt(1.1*2.1) ...
        sqrt(2.1*3.3) sqrt(3.3*4.7) sqrt(4.7*5.8) sqrt(5.8*9) 12], ...
    'diameter_basis','aerodynamic','mass_basis','dry_UO2F2','source_file','');
end

function values = required_numeric(T,name)
values = optional_numeric(T,name);
if isempty(values)
    error('QES:PSD:Schema','Required aerosol sampler column is missing: %s',name);
end
end

function values = optional_numeric(T,name)
values = getcol(T,{name},false);
if ~isempty(values)
    if ~isnumeric(values) && ~islogical(values)
        values = str2double(to_cellstr_compat(values));
    end
    values = double(values(:));
end
end

function values = required_text(T,name)
values = getcol(T,{name},false);
if isempty(values)
    error('QES:PSD:Schema','Required aerosol sampler column is missing: %s',name);
end
values = to_cellstr_compat(values);
values = values(:);
end

function validate_fraction(values)
total = sum(values);
if total ~= 0 && abs(total-1) > 1e-8
    error('QES:PSD:Schema','A nonempty distribution must sum to one.');
end
end

function consistent_statistic(values)
if all(isnan(values)), return; end
if any(~isfinite(values)) || any(values <= 0) ...
        || any(abs(values-values(1)) > 1e-10*max(1,abs(values(1))))
    error('QES:PSD:Schema','Inconsistent PSD summary statistic.');
end
end
