function out = check_CESTA1987_inputs()
cfg = CESTA1987_config();
out = struct();

fprintf('Source XYZ [m]         = %.3f %.3f %.3f\n',cfg.sourceXYZ);
fprintf('Source mass [kg]       = %.3f\n',cfg.sourceMass_kg);
fprintf('Release window [s]     = %.1f -> %.1f (%.1f s)\n',...
    cfg.releaseStart_s,cfg.releaseEnd_s,cfg.releaseDuration_s);
fprintf('Release rate [g/s]     = %.3f\n',cfg.releaseRate_g_s);
fprintf('Exit velocity [m/s]    = %.3f\n',cfg.exitVelocity_m_s);
fprintf('Ambient T/RH/P         = %.2f K / %.3f / %.0f Pa\n',...
    cfg.temperature_K,cfg.relativeHumidity,cfg.pressure_Pa);
fprintf('Stability              = %s\n',cfg.stabilityClass);
fprintf('Wind profile flag      = direct data-entry profile (flag 4)\n');
fprintf('Wind z [m]             = %s\n',mat2str(cfg.windZ_m'));
fprintf('Wind speed [m/s]       = %s\n',mat2str(cfg.windSpeed_m_s'));
fprintf('QES wind-from azimuth  = %.1f deg\n',cfg.windFrom_deg);
fprintf('CESTA plume azimuth    = %.1f deg\n',cfg.plumeAzimuth_deg);

names = {'observation_sampler.csv','component_history.csv',...
         'parcel_mass_diagnostic.csv','lmz_windsOut.nc',...
         'lmz_windsWk.nc','lmz_turbOut.nc','lmz_plumeOut.nc'};
existsFlag = false(size(names));
for i=1:numel(names)
    f = fullfile(cfg.outputDir,names{i});
    existsFlag(i) = exist(f,'file') == 2;
    fprintf('%-30s : %s\n',names{i},tf(existsFlag(i)));
end
out.files = table(names(:),existsFlag(:),'VariableNames',{'File','Exists'});
end

function s=tf(x)
if x, s='FOUND'; else, s='MISSING'; end
end
