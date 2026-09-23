function out = analyze_parcel_mass_CESTA1987()
cfg = CESTA1987_config();
out = struct();
if ~exist(cfg.parcelCsv,'file')
    fprintf('Parcel diagnostic not found: %s\n',cfg.parcelCsv);
    return;
end

T = readtable_compat(cfg.parcelCsv);
pid = getcol(T,{'particleID','particleId'});
time = getcol(T,{'time'});

qold = getcol(T,{'qMass_old','qMassOld'},false);
qnew = getcol(T,{'qMass_new','qMassNew'},false);
E    = getcol(T,{'E'},false);
us   = getcol(T,{'us','usold'},false);
dt   = getcol(T,{'par_dt','parDt'},false);

fprintf('Parcel diagnostic: rows=%d trackedParticles=%d\n',size(T,1),numel(unique(pid)));

if ~isempty(qold) && ~isempty(qnew) && ~isempty(E) && ~isempty(us) && ~isempty(dt)
    r = qnew-qold-E.*abs(us).*dt;
    den=max(abs(qnew),1e-30);
    rr=abs(r)./den;
    out.maxAbsQMassResidual=max(abs(r(isfinite(r))));
    out.maxRelativeQMassResidual=max(rr(isfinite(rr)));
    out.countRelAbove1e12=sum(rr>1e-12 & isfinite(rr));
    fprintf('qMass accepted-step reconstruction:\n');
    fprintf('  max absolute residual = %.6e\n',out.maxAbsQMassResidual);
    fprintf('  max relative residual = %.6e\n',out.maxRelativeQMassResidual);
    fprintf('  count rel > 1e-12     = %d\n',out.countRelAbove1e12);

    fig=figure('Name','qMass residual','Color','w');
    semilogy(time,max(rr,1e-30),'.','MarkerSize',4);
    grid on; xlabel('Model time (s)'); ylabel('Relative qMass residual');
    title('Accepted-step qMass reconstruction residual');
    savefig_png(fig,fullfile(cfg.resultsDir,'parcel_qMass_relative_residual.png'));
else
    fprintf('qMass_old/qMass_new/E/us/par_dt columns not all found; basic inventory only.\n');
end
end
