function results = run_all_CESTA1987()
% Run the complete CESTA 1987 post-processing suite.
clc;
cfg = CESTA1987_config();

fprintf('\n============================================================\n');
fprintf('QES-Dense CESTA 1987 analysis suite\n');
fprintf('Project : %s\n', cfg.projectDir);
fprintf('Output  : %s\n', cfg.outputDir);
fprintf('Results : %s\n', cfg.resultsDir);
fprintf('============================================================\n\n');

results = struct();

fprintf('[1/9] Input/output consistency check...\n');
results.inputs = check_CESTA1987_inputs();

fprintf('\n[2/9] Wind-profile / wind-NetCDF analysis...\n');
results.wind = analyze_wind_CESTA1987();

fprintf('\n[3/9] Solved source-point wind-profile verification...\n');
results.windSolved = check_solved_wind_profile_CESTA1987();

fprintf('\n[4/9] Component/reaction history analysis...\n');
results.component = analyze_component_history_CESTA1987();

fprintf('\n[5/9] Parcel-mass diagnostic analysis...\n');
results.parcel = analyze_parcel_mass_CESTA1987();

fprintf('\n[6/9] Integrated sampler comparison...\n');
results.integrated = analyze_integrated_CESTA1987();

fprintf('\n[7/9] Sequential sampler comparison...\n');
results.sequential = analyze_sequential_CESTA1987();

fprintf('\n[8/9] Oak Ridge/QES/CESTA PSD evolution analysis...\n');
results.psd = analyze_aerosol_PSD_CESTA1987();

fprintf('\n[9/9] Optional plume NetCDF inventory...\n');
results.plumeNetCDF = analyze_plume_netcdf_CESTA1987();

fprintf('\n============================================================\n');
fprintf('CESTA 1987 analysis finished.\n');
fprintf('See: %s\n', cfg.resultsDir);
fprintf('============================================================\n');
end
