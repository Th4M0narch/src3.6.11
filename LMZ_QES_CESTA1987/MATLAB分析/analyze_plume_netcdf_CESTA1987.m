function out = analyze_plume_netcdf_CESTA1987()
% Optional NetCDF inventory.
%
% IMPORTANT:
% The current D2 design keeps UF6/HF/UO2F2 species in the sparse observation
% sampler CSV and does NOT expand the frozen plume NetCDF schema. Therefore
% a legacy NetCDF variable named "conc" must not automatically be interpreted
% as post-reaction UF6/HF/UO2F2 concentration.
cfg=CESTA1987_config();
out=struct();
if ~exist(cfg.plumeNc,'file')
    fprintf('Plume NetCDF not found: %s\n',cfg.plumeNc);
    return;
end
info=ncinfo(cfg.plumeNc);
f=fullfile(cfg.resultsDir,'plume_netcdf_inventory.txt');
fid=fopen(f,'w');
fprintf(fid,'File: %s\n\n',cfg.plumeNc);
for i=1:numel(info.Dimensions)
    fprintf(fid,'DIM %-28s %d\n',info.Dimensions(i).Name,info.Dimensions(i).Length);
end
fprintf(fid,'\n');
for i=1:numel(info.Variables)
    fprintf(fid,'VAR %-28s size=%s\n',info.Variables(i).Name,mat2str(info.Variables(i).Size));
end
fclose(fid);
fprintf('Plume NetCDF inventory written: %s\n',f);
out.variableNames={info.Variables.Name};
end
