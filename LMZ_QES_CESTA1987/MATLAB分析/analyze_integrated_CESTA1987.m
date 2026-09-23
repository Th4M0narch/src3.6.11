function out = analyze_integrated_CESTA1987()
cfg = CESTA1987_config();
out = struct();
if ~exist(cfg.observationCsv,'file')
    fprintf('Observation sampler output not found: %s\n',cfg.observationCsv);
    return;
end

T=readtable_compat(cfg.observationCsv);
typ=to_cellstr_compat(getcol(T,{'sampler_type','type'}));
q=strcmpi(typ,'integrated');
Ti=T(q,:);
if isempty(Ti), fprintf('No integrated rows found.\n'); return; end

dist=getcol(Ti,{'distance'});
azi =getcol(Ti,{'azimuth'});
Up  =getcol(Ti,{'C_U_pred_mg_m3'});
Fp  =getcol(Ti,{'C_F_pred_mg_m3'});
Um  =getcol(Ti,{'C_U_meas_mg_m3','measuredU'});
Fm  =getcol(Ti,{'C_F_meas_mg_m3','measuredF'});
Rp  =getcol(Ti,{'R_UF_pred'},false);
Rm  =getcol(Ti,{'R_UF_meas'},false);

MU=cesta_metrics(Up,Um);
MF=cesta_metrics(Fp,Fm);
fprintf('Integrated U: N=%d FAC2=%.3f FAC5=%.3f FAC10=%.3f MG=%.3f VG=%.3f FB=%.3f NMSE=%.3f NAD=%.3f CC=%.3f\n',...
    MU.N,MU.FAC2,MU.FAC5,MU.FAC10,MU.MG,MU.VG,MU.FB,MU.NMSE,MU.NAD,MU.CC);
fprintf('Integrated F: N=%d FAC2=%.3f FAC5=%.3f FAC10=%.3f MG=%.3f VG=%.3f FB=%.3f NMSE=%.3f NAD=%.3f CC=%.3f\n',...
    MF.N,MF.FAC2,MF.FAC5,MF.FAC10,MF.MG,MF.VG,MF.FB,MF.NMSE,MF.NAD,MF.CC);

MT=[metrics_to_table('Integrated_U',MU); metrics_to_table('Integrated_F',MF)];
writetable(MT,fullfile(cfg.resultsDir,'integrated_metrics.csv'));

scatter_factor_plot(Um,Up,'Integrated uranium: model vs CESTA',...
    'Measured U (mg/m^3)','Predicted U (mg/m^3)',...
    fullfile(cfg.resultsDir,'integrated_U_scatter.png'));
scatter_factor_plot(Fm,Fp,'Integrated fluorine: model vs CESTA',...
    'Measured F (mg/m^3)','Predicted F (mg/m^3)',...
    fullfile(cfg.resultsDir,'integrated_F_scatter.png'));

% Crosswind profiles by radius.
radii=unique(dist(isfinite(dist)))';
for r=radii
    k=dist==r & isfinite(azi);
    [az,ix]=sort(azi(k));
    up=Up(k); up=up(ix); um=Um(k); um=um(ix);
    fp=Fp(k); fp=fp(ix); fm=Fm(k); fm=fm(ix);

    fig=figure('Color','w');
    subplot(2,1,1);
    semilogy(az,um,'ko-',az,up,'o--','LineWidth',1.0);
    grid on; xlabel('Azimuth (deg)'); ylabel('U (mg/m^3)');
    title(sprintf('Integrated U, r=%g m',r)); legend('Measured','Predicted','Location','best');
    subplot(2,1,2);
    semilogy(az,fm,'ko-',az,fp,'o--','LineWidth',1.0);
    grid on; xlabel('Azimuth (deg)'); ylabel('F (mg/m^3)');
    title(sprintf('Integrated F, r=%g m',r)); legend('Measured','Predicted','Location','best');
    savefig_png(fig,fullfile(cfg.resultsDir,sprintf('integrated_crosswind_R%g.png',r)));
end

% Downwind maxima.
R=[]; UmaxM=[]; UmaxP=[]; FmaxM=[]; FmaxP=[];
for r=radii
    k=dist==r;
    R(end+1,1)=r;
    UmaxM(end+1,1)=safeMax(Um(k));
    UmaxP(end+1,1)=safeMax(Up(k));
    FmaxM(end+1,1)=safeMax(Fm(k));
    FmaxP(end+1,1)=safeMax(Fp(k));
end
MaxT=table(R,UmaxM,UmaxP,FmaxM,FmaxP,...
    'VariableNames',{'distance_m','U_meas_max','U_pred_max','F_meas_max','F_pred_max'});
writetable(MaxT,fullfile(cfg.resultsDir,'integrated_downwind_maxima.csv'));

fig=figure('Color','w');
loglog(R,UmaxM,'ko-',R,UmaxP,'o--',R,FmaxM,'ks-',R,FmaxP,'s--','LineWidth',1.0);
grid on; xlabel('Downwind distance (m)'); ylabel('Maximum concentration (mg/m^3)');
title('CESTA integrated-sampler downwind maxima');
legend('U measured','U predicted','F measured','F predicted','Location','best');
savefig_png(fig,fullfile(cfg.resultsDir,'integrated_downwind_maxima.png'));

% Normalized U/F ratio.
if isempty(Rp)
    Rp=(114/238).*Up./Fp;
end
if isempty(Rm)
    Rm=(114/238).*Um./Fm;
end
qr=isfinite(Rp)&isfinite(Rm)&Rp>=0&Rm>=0;
if any(qr)
    fig=figure('Color','w');
    plot(Rm(qr),Rp(qr),'o'); hold on;
    lo=min([Rm(qr);Rp(qr)]); hi=max([Rm(qr);Rp(qr)]);
    plot([lo hi],[lo hi],'k-'); grid on;
    xlabel('Measured normalized U/F ratio'); ylabel('Predicted normalized U/F ratio');
    title('Integrated normalized U/F ratio');
    savefig_png(fig,fullfile(cfg.resultsDir,'integrated_RUF.png'));
end

out.U=MU; out.F=MF; out.maxima=MaxT; out.rows=height(Ti);
end

function m=safeMax(x)
x=x(isfinite(x));
if isempty(x), m=NaN; else, m=max(x); end
end
