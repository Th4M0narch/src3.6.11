function plot_aerosol_process_diagnostics(D,outputDir)
% Plot deposition, coagulation, mass-balance, and distance diagnostics.

if nargin < 2 || isempty(outputDir)
    cfg = CESTA1987_config();
    outputDir = fullfile(cfg.outputDir,'PSD_analysis');
end
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

plot_deposition_fraction(D.distance,outputDir);
plot_coagulation_growth(D.distance,outputDir);
plot_distance_growth_loss(D.distance,outputDir);
plot_mass_balance(D.massBalance,outputDir);
end

function plot_deposition_fraction(T,outputDir)
fig = figure('Color','w');
hold on;
hasData = false;
if ~isempty(T) && size(T,1) > 0
    q = isfinite(T.deposition_fraction);
    if any(q)
        plot(mean_distance(T(q,:)),100.0*T.deposition_fraction(q), ...
            'ko-','LineWidth',1.2,'MarkerSize',5);
        hasData = true;
    end
end
grid on;
xlabel('Distance (m)');
ylabel('UO2F2 deposition fraction (%)');
title('Aerosol deposition fraction by distance');
if hasData
    legend('QES deposited inventory','Location','best');
else
    text(0.1,0.5,'No deposition data','Units','normalized');
end
savefig_png(fig,fullfile(outputDir,'deposition_fraction_distance.png'));
close(fig);
end

function plot_coagulation_growth(T,outputDir)
fig = figure('Color','w');
hold on;
hasData = false;
if ~isempty(T) && size(T,1) > 0
    q = isfinite(T.coagulation_growth);
    if any(q)
        plot(mean_distance(T(q,:)),T.coagulation_growth(q), ...
            'bs-','LineWidth',1.2,'MarkerSize',5);
        hasData = true;
    end
end
grid on;
xlabel('Distance (m)');
ylabel('Coagulation diameter growth factor');
title('Coagulation growth by distance');
if hasData
    legend('Mean cumulative growth','Location','best');
else
    text(0.1,0.5,'No coagulation data','Units','normalized');
end
savefig_png(fig,fullfile(outputDir,'coagulation_growth_distance.png'));
close(fig);
end

function plot_distance_growth_loss(T,outputDir)
fig = figure('Color','w');
hold on;
hasData = false;
if ~isempty(T) && size(T,1) > 0
    x = mean_distance(T);
    if any(isfinite(T.d50_growth))
        plot(x,T.d50_growth,'bs-','LineWidth',1.2,'MarkerSize',5);
        hasData = true;
    end
    if any(isfinite(T.mass_loss_fraction))
        plot(x,T.mass_loss_fraction,'ko--','LineWidth',1.2,'MarkerSize',5);
        hasData = true;
    end
end
grid on;
xlabel('Distance (m)');
ylabel('Ratio');
title('Distance-direction d50 growth and aerosol mass loss');
if hasData
    legend('d50 growth factor','Mass loss fraction','Location','best');
else
    text(0.1,0.5,'No distance diagnostics','Units','normalized');
end
savefig_png(fig,fullfile(outputDir,'d50_growth_mass_loss_distance.png'));
close(fig);
end

function plot_mass_balance(T,outputDir)
fig = figure('Color','w');
hold on;
hasData = false;
if ~isempty(T) && size(T,1) > 0
    q = isfinite(T.time_s);
    if any(q)
        plot(T.time_s(q),T.uo2f2_before_kg(q),'k-','LineWidth',1.1);
        plot(T.time_s(q),T.uo2f2_after_deposition_kg(q),'b-','LineWidth',1.1);
        plot(T.time_s(q),T.uo2f2_deposited_kg(q),'r--','LineWidth',1.1);
        hasData = true;
    end
end
grid on;
xlabel('Time (s)');
ylabel('UO2F2 mass (kg)');
title('Aerosol mass balance by output time');
if hasData
    legend('Before aerosol step','Airborne after deposition', ...
        'Deposited during step','Location','best');
else
    text(0.1,0.5,'No mass-balance data','Units','normalized');
end
savefig_png(fig,fullfile(outputDir,'aerosol_mass_balance_time.png'));
close(fig);
end

function x = mean_distance(T)
x = 0.5 * (T.distance_lower_m + T.distance_upper_m);
q = isfinite(x);
if ~all(q)
    x(~q) = T.distance_lower_m(~q);
end
end
