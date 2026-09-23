function out = analyze_aerosol_PSD_CESTA1987()
% Analyze QES aerosol PSD evolution and compare it with the Oak Ridge PSD.
%
% The QES sampler CSV stores window-accumulated mass/number fractions for
% each sampler and output time.  This script preserves the file rows and
% reports both the raw QES evolution and metrics against the experimental
% reference PSD.  Missing reference data does not prevent QES-only output.

cfg = CESTA1987_config();
out = struct();
out.status = 'not_started';
out.outputDir = fullfile(cfg.outputDir,'PSD_analysis');
out.metricsFile = fullfile(out.outputDir,'psd_metrics.csv');
out.binComparisonFile = fullfile(out.outputDir,'psd_bin_comparison.csv');
if ~exist(out.outputDir,'dir')
    mkdir(out.outputDir);
end

fprintf('Reading Oak Ridge initial aerosol PSD...\n');
psd = read_oakridge_psd_analysis();
if isempty(psd.exp.diameter_um)
    warning(['Oak Ridge reference PSD unavailable. ' ...
        'QES PSD evolution will be produced without experimental metrics.']);
end

fprintf('Reading QES aerosol PSD sampler output...\n');
qes = read_qes_psd_sampler_analysis();
psd.qes = qes.qes;

if isempty(psd.qes)
    warning('No QES aerosol PSD sampler rows found. PSD analysis skipped.');
    out.process = analyze_aerosol_process_diagnostics(out.outputDir);
    out.status = 'missing_qes_psd';
    out.psd = psd;
    return;
end

cfg = CESTA1987_config();
metadata = load_cesta_sampler_metadata(cfg);
matchingReport = table();
if ~isempty(metadata)
    [~,matchingReport] = match_cesta_sampler(metadata,psd.qes);
    writetable(matchingReport,fullfile(out.outputDir,'matching_report.csv'));
else
    warning(['CESTA sampler metadata unavailable; QES PSD evolution will ' ...
        'be reported without experimental sampler matches.']);
end

% Oak Ridge remains the global initial condition reference.  It must not
% be copied to every CESTA sampler as if it were a local experimental PSD.
% With no per-sampler experimental PSD table, comparison metrics stay NaN.
metrics = calculate_aerosol_psd_metrics(psd,struct());
metricsTable = metrics_to_output_table(metrics);
writetable(metricsTable,out.metricsFile);

binTable = bin_comparison_table(psd,matchingReport);
writetable(binTable,out.binComparisonFile);

plot_aerosol_psd_results(psd,metrics,out.outputDir);
out.process = analyze_aerosol_process_diagnostics(out.outputDir);

fprintf('Aerosol PSD analysis complete.\n');
fprintf('  QES PSD rows: %d\n',numel(psd.qes));
fprintf('  Oak Ridge reference: %s\n',ternary(isempty(psd.exp.diameter_um), ...
    'unavailable','available'));
fprintf('  Output: %s\n',out.outputDir);

out.status = 'complete';
out.psd = psd;
out.metrics = metricsTable;
out.binComparison = binTable;
out.matchingReport = matchingReport;
end

function T = metrics_to_output_table(metrics)
T = table( ...
    metrics.distance_m, ...
    metrics.height_m, ...
    metrics.time_s, ...
    metrics.d50_um, ...
    metrics.sigma_g, ...
    metrics.d50_change, ...
    metrics.sigma_change, ...
    metrics.mass_RMSE, ...
    metrics.number_RMSE, ...
    metrics.JS_divergence, ...
    'VariableNames', ...
    {'distance_m','height_m','time_s','d50_um','sigma_g', ...
     'd50_change','sigma_change','mass_RMSE','number_RMSE', ...
     'JS_divergence'});
end

function T = bin_comparison_table(psd,matchingReport)
n = numel(psd.qes);
rows = n * 9;
distance_m = NaN(rows,1);
height_m = NaN(rows,1);
time_s = NaN(rows,1);
bin = NaN(rows,1);
reference_mass_fraction = NaN(rows,1);
qes_mass_fraction = NaN(rows,1);
difference = NaN(rows,1);

row = 0;
for i = 1:n
    for b = 1:9
        row = row + 1;
        distance_m(row) = psd.qes(i).distance;
        height_m(row) = psd.qes(i).height;
        time_s(row) = psd.qes(i).time;
        bin(row) = b;
        qes_mass_fraction(row) = psd.qes(i).mass_fraction(b);
        reference_mass_fraction(row) = NaN;
        difference(row) = NaN;
    end
end

T = table(distance_m,height_m,time_s,bin, ...
    reference_mass_fraction,qes_mass_fraction,difference);
if ~isempty(matchingReport)
    T.Properties.UserData = matchingReport;
end
end

function metadata = load_cesta_sampler_metadata(cfg)
metadata = struct('id',{},'type',{},'distance',{},'height',{}, ...
    'azimuth',{},'windowStart',{},'windowEnd',{});
if ~exist(cfg.observationCsv,'file')
    return;
end

T = readtable_compat(cfg.observationCsv);
n = size(T,1);
id = to_cellstr_compat(getcol(T,{'sampler_id','id'},false));
typ = to_cellstr_compat(getcol(T,{'sampler_type','type'},false));
dist = getcol(T,{'distance','distance_m'},false);
height = getcol(T,{'z','height','height_m'},false);
azi = getcol(T,{'azimuth','azimuth_deg'},false);
w0 = getcol(T,{'window_start','windowStart'},false);
w1 = getcol(T,{'window_end','windowEnd'},false);
if isempty(dist) || isempty(height)
    return;
end

for k = 1:n
    if isempty(id)
        idText = sprintf('sampler_%d',k);
    else
        idText = id{k};
    end
    if isempty(typ)
        typeText = 'sampler';
    else
        typeText = typ{k};
    end
    azimuth = NaN;
    if ~isempty(azi)
        azimuth = azi(k);
    end
    windowStart = cfg.releaseStart_s;
    windowEnd = cfg.releaseEnd_s;
    if ~isempty(w0), windowStart = w0(k); end
    if ~isempty(w1), windowEnd = w1(k); end
    entry = struct('id',char(idText),'type',char(typeText), ...
        'distance',double(dist(k)),'height',double(height(k)), ...
        'azimuth',double(azimuth),'windowStart',double(windowStart), ...
        'windowEnd',double(windowEnd));
    metadata(end+1) = entry; %#ok<AGROW>
end
end

function value = ternary(condition,trueValue,falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
