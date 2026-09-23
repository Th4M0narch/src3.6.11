function out = analyze_component_history_CESTA1987()
cfg = CESTA1987_config();
out = struct();
if ~exist(cfg.componentCsv,'file')
    fprintf('Component history not found: %s\n',cfg.componentCsv);
    return;
end

T = readtable_compat(cfg.componentCsv);

pid  = getcol(T,{'particleID','particleId'});
age  = getcol(T,{'particleAge'},false);
time = getcol(T,{'time'});
relT = getcol(T,{'releaseTime'},false);
if isempty(age)
    if isempty(relT), error('Need particleAge or releaseTime in component_history.csv'); end
    age = time-relT;
end

mAir = getcol(T,{'mAirCum'});
mP   = getcol(T,{'mPrimary'});
mG   = getcol(T,{'mGasProduct'});
mPart= getcol(T,{'mParticleProduct'});
mW   = getcol(T,{'mAmbientWater'});
conv = getcol(T,{'conversionFraction'});
Rmass= getcol(T,{'Rmass'});
relRU= getcol(T,{'relRU'});
relRF= getcol(T,{'relRF'});

ids = unique(pid(:))';
fprintf('Component history: rows=%d trackedParticles=%d\n',size(T,1),numel(ids));
fprintf('max |Rmass| = %.6e kg\n',max(abs(Rmass(isfinite(Rmass)))));
fprintf('max |relRU| = %.6e\n',max(abs(relRU(isfinite(relRU)))));
fprintf('max |relRF| = %.6e\n',max(abs(relRF(isfinite(relRF)))));

fig=figure('Name','C3 conversion history','Color','w');
hold on;
for k=1:numel(ids)
    q=pid==ids(k);
    [aa,ix]=sort(age(q));
    cc=conv(q); cc=cc(ix);
    plot(aa,cc,'LineWidth',0.8);
end
grid on; xlabel('Particle age (s)'); ylabel('Conversion fraction');
title('CESTA 1987 tracked-particle conversion histories');
savefig_png(fig,fullfile(cfg.resultsDir,'component_conversion_history.png'));

fig=figure('Name','C3 entrained air','Color','w');
hold on;
for k=1:numel(ids)
    q=pid==ids(k);
    [aa,ix]=sort(age(q));
    mm=mAir(q); mm=mm(ix);
    plot(aa,mm,'LineWidth',0.8);
end
grid on; xlabel('Particle age (s)'); ylabel('mAirCum (kg)');
title('Cumulative entrained-air mass');
savefig_png(fig,fullfile(cfg.resultsDir,'component_mAirCum_history.png'));

% Fractions relative to each particle's initial primary mass.
fig=figure('Name','C3 component fractions','Color','w');
hold on;
summary = nan(numel(ids),8);
for k=1:numel(ids)
    q=find(pid==ids(k));
    [aa,ix]=sort(age(q)); q=q(ix);
    rr=getcol(T(q,:),{'reactedPrimary'},false);
    if isempty(rr)
        m0=max(mP(q));
    else
        rr=rr(isfinite(rr));
        % Complete-history reconstruction: final remaining primary
        % plus all accepted-step reacted primary.
        m0=mP(q(end)) + sum(rr);
    end
    if ~isfinite(m0) || m0<=0, m0=max(mP(q)); end
    plot(aa,mP(q)/m0,'LineWidth',0.7);
    summary(k,:)=[ids(k),max(aa),m0,mP(q(end)),mG(q(end)),mPart(q(end)),mW(q(end)),conv(q(end))];
end
grid on; xlabel('Particle age (s)'); ylabel('mPrimary/mPrimary0');
title('Remaining primary fraction, 20 tracked particles');
savefig_png(fig,fullfile(cfg.resultsDir,'component_primary_fraction.png'));

S = array2table(summary,'VariableNames',...
    {'particleID','maxAge_s','mPrimary0_kg','mPrimaryFinal_kg',...
     'mGasFinal_kg','mParticleFinal_kg','mWaterFinal_kg','conversionFinal'});
writetable(S,fullfile(cfg.resultsDir,'component_particle_summary.csv'));

out.rows = size(T,1);
out.trackedParticles = numel(ids);
out.maxAbsRmass = max(abs(Rmass(isfinite(Rmass))));
out.maxAbsRelRU = max(abs(relRU(isfinite(relRU))));
out.maxAbsRelRF = max(abs(relRF(isfinite(relRF))));
out.summary = S;
end
