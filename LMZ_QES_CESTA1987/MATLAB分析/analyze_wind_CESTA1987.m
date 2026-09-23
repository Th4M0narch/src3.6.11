function out = analyze_wind_CESTA1987()
cfg = CESTA1987_config();
out = struct();

z = cfg.windZ_m;
U = cfg.windSpeed_m_s;

% Fit a reference power law U = U10*(z/10)^p to the four direct-entry data.
U10 = interp1(z,U,10,'linear');
p = sum(log(z/10).*log(U/U10))/sum(log(z/10).^2);
U_pow = U10*(z/10).^p;

% Fit U = a ln(z) + b, then infer equivalent z0 from b=-a ln(z0).
coef = polyfit(log(z),U,1);
a = coef(1); b = coef(2);
z0eq = exp(-b/a);
U_log = a*log(z/z0eq);

fprintf('Direct-entry wind profile from sensor XML:\n');
for i=1:numel(z)
    fprintf('  z=%5.1f m  U=%6.3f m/s\n',z(i),U(i));
end
fprintf('Reference power-law exponent p = %.6f\n',p);
fprintf('Reference log-fit equivalent z0 = %.6g m\n',z0eq);

fig=figure('Name','CESTA1987 wind profile','Color','w');
plot(U,z,'ko-','LineWidth',1.2); hold on;
plot(U_pow,z,'--','LineWidth',1.1);
plot(U_log,z,':','LineWidth',1.1);
grid on; xlabel('Wind speed (m/s)'); ylabel('Height z (m)');
title('CESTA 1987 direct-entry mean wind profile');
legend('Direct input','Power-law reference','Log-fit reference','Location','best');
savefig_png(fig,fullfile(cfg.resultsDir,'wind_profile_input.png'));

out.powerExponent = p;
out.logEquivalentZ0_m = z0eq;

% NetCDF inventory and optional profile extraction.
ncFile = cfg.windsVizNc;
if ~exist(ncFile,'file')
    fprintf('Wind NetCDF not found, skipping field extraction: %s\n',ncFile);
    return;
end

info = ncinfo(ncFile);
invFile = fullfile(cfg.resultsDir,'wind_netcdf_inventory.txt');
fid=fopen(invFile,'w');
fprintf(fid,'File: %s\n\n',ncFile);
for i=1:numel(info.Variables)
    fprintf(fid,'%-32s size=%s\n',info.Variables(i).Name,mat2str(info.Variables(i).Size));
end
fclose(fid);
fprintf('NetCDF inventory written: %s\n',invFile);

names = {info.Variables.Name};
uName = findName(names,{'u','U','uMean','u_velocity','uVel'});
vName = findName(names,{'v','V','vMean','v_velocity','vVel'});
zName = findName(names,{'z','z_cc','zCell','z_center','z_face'});

if isempty(uName) || isempty(vName)
    fprintf('Could not unambiguously identify collocated u/v variables; inventory only.\n');
    return;
end

try
    uu = double(ncread(ncFile,uName));
    vv = double(ncread(ncFile,vName));
    if isequal(size(uu),size(vv))
        spd = sqrt(uu.^2+vv.^2);
        out.windSpeedMin = min(spd(:));
        out.windSpeedMax = max(spd(:));
        out.windSpeedMean = mean(spd(isfinite(spd)));
        fprintf('Wind field |U_h|: mean=%.4f min=%.4f max=%.4f m/s\n',...
            out.windSpeedMean,out.windSpeedMin,out.windSpeedMax);
    else
        fprintf('u/v are staggered with different dimensions; no direct magnitude calculation.\n');
    end
catch ME
    fprintf('Wind NetCDF extraction skipped: %s\n',ME.message);
end
end

function n=findName(names,cands)
n='';
for i=1:numel(cands)
    k=find(strcmpi(names,cands{i}),1);
    if ~isempty(k), n=names{k}; return; end
end
end
