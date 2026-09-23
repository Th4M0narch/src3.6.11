function out = analyze_aerosol_process_diagnostics(outputDir)
% Analyze aerosol deposition, coagulation growth, and mass closure.

if nargin < 1 || isempty(outputDir)
    cfg = CESTA1987_config();
    outputDir = fullfile(cfg.outputDir,'PSD_analysis');
end
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

out = struct();
out.status = 'not_started';
out.outputDir = outputDir;
inputs = read_aerosol_process_inputs();
if ~inputs.hasMass
    warning('Aerosol process diagnostics skipped: mass diagnostic missing.');
    out.status = 'missing_mass_diagnostic';
    out.inputs = inputs;
    return;
end

D = calculate_aerosol_process_diagnostics(inputs);
if isempty(D.particle)
    warning('Aerosol process diagnostics contain no particle rows.');
    out.status = 'no_particle_diagnostics';
    out.inputs = inputs;
    return;
end

out.particleFile = fullfile(outputDir,'psd_process_diagnostics.csv');
out.distanceFile = fullfile(outputDir,'psd_distance_diagnostics.csv');
out.coagulationFile = fullfile(outputDir,'coagulation_growth_summary.csv');
out.massBalanceFile = fullfile(outputDir,'aerosol_mass_balance.csv');
out.massSummaryFile = fullfile(outputDir,'aerosol_mass_balance_summary.csv');
out.depositionFile = fullfile(outputDir,'aerosol_deposition_bin_summary.csv');

writetable(D.particle,out.particleFile);
writetable(D.distance,out.distanceFile);
writetable(D.coagulation,out.coagulationFile);
writetable(D.massBalance,out.massBalanceFile);
writetable(D.massBalanceSummary,out.massSummaryFile);
if ~isempty(D.depositionBins)
    writetable(D.depositionBins,out.depositionFile);
end

plot_aerosol_process_diagnostics(D,outputDir);

fprintf('Aerosol process diagnostics complete.\n');
fprintf('  Tracked particles: %d\n',height(D.particle));
fprintf('  Mass-balance rows: %d\n',height(D.massBalance));
fprintf('  Distance groups: %d\n',height(D.distance));

out.status = 'complete';
out.inputs = inputs;
out.results = D;
end
