function out = analyze_psd_CESTA1987()
% Compare the Oak Ridge initial PSD with QES PSD output at CESTA samplers.

cfg = CESTA1987_config();
out = struct();
out.status = 'not_started';
out.outputDir = fullfile(cfg.outputDir,'psd_analysis');
if ~exist(out.outputDir,'dir')
    mkdir(out.outputDir);
end

fprintf('Reading Oak Ridge initial PSD...\n');
psd = struct();
psd.exp = read_oakridge_psd();
if isempty(psd.exp.diameter_um)
    warning('Oak Ridge PSD unavailable; QES-only PSD analysis will continue.');
end

fprintf('Reading QES PSD sampler output...\n');
psd.qes = read_qes_psd_sampler();
if isempty(psd.qes.time_s)
    warning('No QES PSD sampler rows found. PSD analysis skipped.');
    out.status = 'missing_qes_psd';
    return;
end

metadata = load_sampler_metadata(cfg,psd.qes);
if isempty(metadata)
    warning('No CESTA sampler metadata found. PSD analysis skipped.');
    out.status = 'missing_sampler_metadata';
    return;
end

[matchTable,matchingReport] = match_cesta_sampler(metadata,psd.qes);
writetable(matchingReport,fullfile(out.outputDir,'matching_report.csv'));

n = height(matchTable);
qesRow = matchTable.qes_row;
validMatch = isfinite(qesRow);
matchDistance = NaN(n,1);
matchHeight = NaN(n,1);
matchAzimuth = NaN(n,1);
matchTime = NaN(n,1);
matchD50 = NaN(n,1);
matchSigma = NaN(n,1);
matchedMass = NaN(n,9);
matchedNumber = NaN(n,9);
for k = 1:n
    if ~validMatch(k)
        continue;
    end
    row = qesRow(k);
    matchDistance(k) = psd.qes.distance_m(row);
    matchHeight(k) = psd.qes.height_m(row);
    matchAzimuth(k) = psd.qes.azimuth_deg(row);
    matchTime(k) = psd.qes.time_s(row);
    matchD50(k) = psd.qes.d50_um(row);
    matchSigma(k) = psd.qes.sigma_g(row);
    matchedMass(k,:) = psd.qes.mass_fraction(row,:);
    matchedNumber(k,:) = psd.qes.number_fraction(row,:);
end

if ~any(validMatch)
    warning('No QES PSD rows could be matched to CESTA sampler positions.');
    out.status = 'no_matches';
    return;
end

% Oak Ridge is only an initial-condition reference.  It is not a measured
% PSD at every CESTA sampler.  Until per-sampler experimental PSD data are
% supplied, experimental comparison fields remain NaN.
expD50 = NaN(n,1);
expSigma = NaN(n,1);
expMass = NaN(n,9);
expNumber = NaN(n,9);

metrics = calculate_psd_metrics( ...
    expMass,matchedMass,expNumber,matchedNumber, ...
    expD50,matchD50,expSigma,matchSigma);

matchedT = table( ...
    matchTable.experiment_sampler_id,matchTable.qes_sampler_id, ...
    matchDistance,matchHeight,matchAzimuth,matchTime, ...
    expD50,matchD50,expSigma,matchSigma, ...
    'VariableNames', ...
    {'experiment_sampler_id','qes_sampler_id','distance','height','azimuth','time', ...
     'd50_exp','d50_qes','sigma_exp','sigma_qes'});

metricsT = table( ...
    matchDistance,matchHeight,matchAzimuth,matchTime, ...
    expD50,matchD50,metrics.d50_error_relative, ...
    expSigma,matchSigma,metrics.sigma_error, ...
    metrics.mass_RMSE,metrics.number_RMSE,metrics.JS_divergence, ...
    'VariableNames', ...
    {'distance','height','azimuth','time', ...
     'd50_exp','d50_qes','d50_error', ...
     'sigma_exp','sigma_qes','sigma_error', ...
     'mass_RMSE','number_RMSE','JS_divergence'});
matchedT = matchedT(validMatch,:);
metricsT = metricsT(validMatch,:);
writetable(matchedT,fullfile(out.outputDir,'matched_psd.csv'));
writetable(metricsT,fullfile(out.outputDir,'psd_metrics.csv'));

nRows = sum(validMatch) * 9;
binDistance = NaN(nRows,1);
binHeight = NaN(nRows,1);
binId = NaN(nRows,1);
expFraction = NaN(nRows,1);
qesFraction = NaN(nRows,1);
difference = NaN(nRows,1);
rowOut = 0;
for k = 1:n
    if ~validMatch(k)
        continue;
    end
    for b = 1:9
        rowOut = rowOut + 1;
        binDistance(rowOut) = matchDistance(k);
        binHeight(rowOut) = matchHeight(k);
        binId(rowOut) = b;
        expFraction(rowOut) = expMass(k,b);
        qesFraction(rowOut) = matchedMass(k,b);
        difference(rowOut) = matchedMass(k,b) - expMass(k,b);
    end
end

binT = table(binDistance,binHeight,binId,expFraction,qesFraction,difference, ...
    'VariableNames', ...
    {'distance','height','bin','exp_mass_fraction','qes_mass_fraction','difference'});
writetable(binT,fullfile(out.outputDir,'psd_bin_comparison.csv'));

plot_psd_comparison(psd,metricsT,out.outputDir);

fprintf('PSD comparison complete.\n');
fprintf('  Samples: %d\n',sum(validMatch));
fprintf('  Mean mass RMSE: %.6g\n',mean_finite(metrics.mass_RMSE));
fprintf('  Mean number RMSE: %.6g\n',mean_finite(metrics.number_RMSE));
fprintf('  Mean JS divergence: %.6g\n',mean_finite(metrics.JS_divergence));
fprintf('  Output: %s\n',out.outputDir);

out.status = 'complete';
out.matched = matchedT;
out.metrics = metricsT;
out.binComparison = binT;
out.psd = psd;
end

function metadata = load_sampler_metadata(cfg,qes)
metadata = empty_metadata();

if isfield(cfg,'integratedSamplers') || isfield(cfg,'sequentialSamplers')
    if isfield(cfg,'integratedSamplers')
        metadata = append_config_samplers(metadata,cfg.integratedSamplers,'integrated');
    end
    if isfield(cfg,'sequentialSamplers')
        metadata = append_config_samplers(metadata,cfg.sequentialSamplers,'sequential');
    end
end

if isempty(metadata) && exist(cfg.observationCsv,'file')
    T = readtable_compat(cfg.observationCsv);
    typ = to_cellstr_compat(getcol(T,{'sampler_type','type'},false));
    dist = getcol(T,{'distance','distance_m'},false);
    height = getcol(T,{'z','height','height_m'},false);
    w0 = getcol(T,{'window_start','windowStart'},false);
    w1 = getcol(T,{'window_end','windowEnd'},false);
    id = to_cellstr_compat(getcol(T,{'sampler_id','id'},false));
    azi = getcol(T,{'azimuth','azimuth_deg'},false);
    if ~isempty(dist) && ~isempty(height)
        for k = 1:numel(dist)
            if isempty(typ)
                typeText = 'sampler';
            else
                typeText = typ{k};
            end
            if isempty(id)
                idText = sprintf('sampler_%d',k);
            else
                idText = id{k};
            end
            if isempty(azi)
                azimuth = NaN;
            else
                azimuth = double(azi(k));
            end
            ws = cfg.releaseStart_s;
            we = cfg.releaseEnd_s;
            if ~isempty(w0), ws = w0(k); end
            if ~isempty(w1), we = w1(k); end
            metadata = append_metadata(metadata,idText,typeText, ...
                double(dist(k)),double(height(k)),azimuth, ...
                double(ws),double(we));
        end
    end
end
end

function metadata = append_config_samplers(metadata,samplers,defaultType)
if isempty(samplers)
    return;
end
if isstruct(samplers)
    for k = 1:numel(samplers)
        s = samplers(k);
        typeText = defaultType;
        idText = sprintf('%s_%d',defaultType,k);
        distance = NaN;
        height = NaN;
        azimuth = NaN;
        windowStart = 0;
        windowEnd = Inf;
        if isfield(s,'type'), typeText = s.type; end
        if isfield(s,'id'), idText = s.id; end
        if isfield(s,'distance'), distance = s.distance; end
        if isfield(s,'height'), height = s.height; end
        if isfield(s,'z'), height = s.z; end
        if isfield(s,'azimuth'), azimuth = s.azimuth; end
        if isfield(s,'azimuth_deg'), azimuth = s.azimuth_deg; end
        if isfield(s,'windowStart'), windowStart = s.windowStart; end
        if isfield(s,'windowEnd'), windowEnd = s.windowEnd; end
        if isfield(s,'window_start'), windowStart = s.window_start; end
        if isfield(s,'window_end'), windowEnd = s.window_end; end
        metadata = append_metadata(metadata,idText,typeText, ...
            distance,height,azimuth,windowStart,windowEnd);
    end
elseif iscell(samplers)
    for k = 1:numel(samplers)
        metadata = append_config_samplers(metadata,samplers{k},defaultType);
    end
end
end

function metadata = append_metadata(metadata,id,type,distance,height,azimuth,w0,w1)
id = char(id);
type = char(type);
entry = struct('id',id,'type',type,'distance',distance,'height',height, ...
    'azimuth',azimuth,'windowStart',w0,'windowEnd',w1);
metadata(end+1) = entry;
end

function metadata = empty_metadata()
metadata = struct('id',{},'type',{},'distance',{},'height',{}, ...
    'azimuth',{},'windowStart',{},'windowEnd',{});
end

function value = mean_finite(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = mean(x);
end
end
