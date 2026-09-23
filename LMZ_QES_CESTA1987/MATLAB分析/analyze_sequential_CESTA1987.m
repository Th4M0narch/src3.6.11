function out = analyze_sequential_CESTA1987()
cfg=CESTA1987_config();
out=struct();
if ~exist(cfg.observationCsv,'file')
    fprintf('Observation sampler output not found: %s\n',cfg.observationCsv);
    return;
end

T=readtable_compat(cfg.observationCsv);
typ=to_cellstr_compat(getcol(T,{'sampler_type','type'}));
q=strcmpi(typ,'sequential');
Ts=T(q,:);
if isempty(Ts), fprintf('No sequential rows found.\n'); return; end

w0=getcol(Ts,{'window_start','windowStart'});
w1=getcol(Ts,{'window_end','windowEnd'});
dur=w1-w0;
dist=getcol(Ts,{'distance'});
azi=getcol(Ts,{'azimuth'});
z=getcol(Ts,{'z'});
Up=getcol(Ts,{'C_U_pred_mg_m3'});
Fp=getcol(Ts,{'C_F_pred_mg_m3'});
Um=getcol(Ts,{'C_U_meas_mg_m3','measuredU'});
Fm=getcol(Ts,{'C_F_meas_mg_m3','measuredF'});

% S1-S7 are about 360 s each; S0 is the 42-min overlapping sampler.
q6 = dur < 500;
MU=cesta_metrics(Up(q6),Um(q6));
MF=cesta_metrics(Fp(q6),Fm(q6));
fprintf('Sequential 6-min U: N=%d FAC2=%.3f FAC5=%.3f FAC10=%.3f MG=%.3f VG=%.3f FB=%.3f NMSE=%.3f NAD=%.3f CC=%.3f\n',...
    MU.N,MU.FAC2,MU.FAC5,MU.FAC10,MU.MG,MU.VG,MU.FB,MU.NMSE,MU.NAD,MU.CC);
fprintf('Sequential 6-min F: N=%d FAC2=%.3f FAC5=%.3f FAC10=%.3f MG=%.3f VG=%.3f FB=%.3f NMSE=%.3f NAD=%.3f CC=%.3f\n',...
    MF.N,MF.FAC2,MF.FAC5,MF.FAC10,MF.MG,MF.VG,MF.FB,MF.NMSE,MF.NAD,MF.CC);

MT=[metrics_to_table('Sequential6min_U',MU); metrics_to_table('Sequential6min_F',MF)];
writetable(MT,fullfile(cfg.resultsDir,'sequential_metrics.csv'));

scatter_factor_plot(Um(q6),Up(q6),'Sequential uranium: model vs CESTA',...
    'Measured U (mg/m^3)','Predicted U (mg/m^3)',...
    fullfile(cfg.resultsDir,'sequential_U_scatter.png'));
scatter_factor_plot(Fm(q6),Fp(q6),'Sequential fluorine: model vs CESTA',...
    'Measured F (mg/m^3)','Predicted F (mg/m^3)',...
    fullfile(cfg.resultsDir,'sequential_F_scatter.png'));

% Time histories for every unique distance/azimuth/height combination.
plotDir=fullfile(cfg.resultsDir,'SequentialPlots');
if ~exist(plotDir,'dir'), mkdir(plotDir); end
key=[dist azi z];
ukey=unique(key,'rows');
for g=1:size(ukey,1)
    r=ukey(g,1); a=ukey(g,2); zz=ukey(g,3);
    k=q6 & dist==r & abs(azi-a)<1e-9 & abs(z-zz)<1e-9;
    if ~any(k), continue; end
    mid=0.5*(w0(k)+w1(k));
    [mid,ix]=sort(mid);
    up=Up(k); up=up(ix); um=Um(k); um=um(ix);
    fp=Fp(k); fp=fp(ix); fm=Fm(k); fm=fm(ix);

    fig=figure('Color','w');
    subplot(2,1,1);
    plot(mid,um,'ko-',mid,up,'o--','LineWidth',1.0);
    grid on; xlabel('Model time (s)'); ylabel('U (mg/m^3)');
    title(sprintf('Sequential U: r=%g m, az=%.1f deg, z=%g m',r,a,zz));
    legend('Measured','Predicted','Location','best');
    subplot(2,1,2);
    plot(mid,fm,'ko-',mid,fp,'o--','LineWidth',1.0);
    grid on; xlabel('Model time (s)'); ylabel('F (mg/m^3)');
    title(sprintf('Sequential F: r=%g m, az=%.1f deg, z=%g m',r,a,zz));
    legend('Measured','Predicted','Location','best');
    fname=sprintf('seq_R%g_A%s_Z%g.png',r,strrep(num2str(a),'.','p'),zz);
    savefig_png(fig,fullfile(plotDir,fname));
    close(fig);
end

% Compare S0 against duration-weighted S1-S7 at the same physical sampler.
S0rows=[];
for g=1:size(ukey,1)
    r=ukey(g,1); a=ukey(g,2); zz=ukey(g,3);
    k0=~q6 & dist==r & abs(azi-a)<1e-9 & abs(z-zz)<1e-9;
    k6=q6  & dist==r & abs(azi-a)<1e-9 & abs(z-zz)<1e-9;
    if sum(k0)==1 && any(k6)
        wu=dur(k6);
        upW=sum(Up(k6).*wu)/sum(wu);
        fpW=sum(Fp(k6).*wu)/sum(wu);
        umW=weightedFinite(Um(k6),wu);
        fmW=weightedFinite(Fm(k6),wu);
        i0=find(k0,1);
        S0rows(end+1,:)=[r a zz Up(i0) upW Fp(i0) fpW Um(i0) umW Fm(i0) fmW]; %#ok<AGROW>
    end
end
if ~isempty(S0rows)
    S0T=array2table(S0rows,'VariableNames',...
       {'distance_m','azimuth_deg','z_m','U_pred_S0','U_pred_weighted_S1S7',...
        'F_pred_S0','F_pred_weighted_S1S7','U_meas_S0','U_meas_weighted_S1S7',...
        'F_meas_S0','F_meas_weighted_S1S7'});
    writetable(S0T,fullfile(cfg.resultsDir,'sequential_S0_vs_S1S7.csv'));
    out.S0Comparison=S0T;
end

out.U=MU; out.F=MF; out.rows=height(Ts); out.rows6min=sum(q6);
end

function y=weightedFinite(x,w)
q=isfinite(x)&isfinite(w)&w>0;
if ~any(q), y=NaN; else, y=sum(x(q).*w(q))/sum(w(q)); end
end
