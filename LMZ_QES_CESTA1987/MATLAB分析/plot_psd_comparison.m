function plot_psd_comparison(psd,metrics,outDir)
% Create the requested Oak Ridge/QES/CESTA PSD comparison figures.

if nargin < 3 || isempty(outDir)
    cfg = CESTA1987_config();
    outDir = fullfile(cfg.outputDir,'psd_analysis');
end
if ~exist(outDir,'dir')
    mkdir(outDir);
end

if ~isfield(psd,'qes') || isempty(psd.qes.time_s)
    warning('No QES PSD rows available for plotting.');
    return;
end

plot_initial(psd,outDir);
plot_distance_evolution(psd,outDir);
plot_d50(metrics,outDir);
plot_sigma(metrics,outDir);
plot_heatmap(psd,outDir);
end

function plot_initial(psd,outDir)
qes = psd.qes;
positive = sum(qes.mass_fraction,2) > 0;
if any(positive)
    t0 = min(qes.time_s(positive));
else
    t0 = min(qes.time_s);
end
q = abs(qes.time_s-t0) <= max(1.0e-12,1.0e-9*abs(t0));
qMass = average_rows(qes.mass_fraction(q,:));
xQes = qes.bin_diameter_um(:);

fig = figure('Color','w');
hold on;
if isfield(psd,'exp') && ~isempty(psd.exp.diameter_um)
    plot(psd.exp.diameter_um(:),psd.exp.mass_fraction(:), ...
        'ko-','LineWidth',1.2,'MarkerSize',5);
end
plot(xQes,qMass,'bs--','LineWidth',1.2,'MarkerSize',5);
set(gca,'XScale','log');
grid on;
xlabel('Geometric diameter (\mum)');
ylabel('Mass fraction per bin');
title(sprintf('Oak Ridge initial PSD vs QES initial snapshot (t = %.3g s)',t0));
if isfield(psd,'exp') && ~isempty(psd.exp.diameter_um)
    legend('Oak Ridge','QES initial','Location','best');
else
    legend('QES initial','Location','best');
end
savefig_png(fig,fullfile(outDir,'oakridge_initial_vs_qes_initial.png'));
close(fig);
end

function plot_distance_evolution(psd,outDir)
targets = [10 20 40 70 100 200 500];
fig = figure('Color','w');
hold on;
legendText = {};
plotted = false;
for k = 1:numel(targets)
    [found,distance_m,distribution] = representative_at_distance(psd.qes,targets(k));
    if found
        plot(psd.qes.bin_diameter_um,distribution,'o-','LineWidth',1.1, ...
            'MarkerSize',4);
        legendText{end+1} = sprintf('%.0f m',distance_m); %#ok<AGROW>
        plotted = true;
    end
end
set(gca,'XScale','log');
grid on;
xlabel('Geometric diameter (\mum)');
ylabel('Mass fraction per bin');
title('QES PSD evolution by nearest CESTA distance');
if plotted
    legend(legendText,'Location','best');
else
    text(0.1,0.5,'No QES PSD rows','Units','normalized');
    axis off;
end
savefig_png(fig,fullfile(outDir,'psd_evolution_distance.png'));
close(fig);
end

function plot_d50(metrics,outDir)
fig = figure('Color','w');
hold on;
hasExp = false;
hasQes = false;
if ~isempty(metrics) && height(metrics) > 0
    [dExp,eExp] = grouped_mean(metrics.distance,metrics.d50_exp);
    [dQes,eQes] = grouped_mean(metrics.distance,metrics.d50_qes);
    if any(isfinite(eExp))
        plot(dExp,eExp,'ko-','LineWidth',1.2,'MarkerSize',5);
        hasExp = true;
    end
    if any(isfinite(eQes))
        plot(dQes,eQes,'bs--','LineWidth',1.2,'MarkerSize',5);
        hasQes = true;
    end
end
grid on;
xlabel('Distance (m)');
ylabel('Mass median diameter d50 (\mum)');
title('d50 evolution');
if hasExp && hasQes
    legend('Experiment (Oak Ridge initial)','QES','Location','best');
elseif hasQes
    legend('QES','Location','best');
elseif hasExp
    legend('Experiment (Oak Ridge initial)','Location','best');
else
    text(0.1,0.5,'No d50 data','Units','normalized');
end
savefig_png(fig,fullfile(outDir,'d50_decay.png'));
close(fig);
end

function plot_sigma(metrics,outDir)
fig = figure('Color','w');
hold on;
hasExp = false;
hasQes = false;
if ~isempty(metrics) && height(metrics) > 0
    [dExp,eExp] = grouped_mean(metrics.distance,metrics.sigma_exp);
    [dQes,eQes] = grouped_mean(metrics.distance,metrics.sigma_qes);
    if any(isfinite(eExp))
        plot(dExp,eExp,'ko-','LineWidth',1.2,'MarkerSize',5);
        hasExp = true;
    end
    if any(isfinite(eQes))
        plot(dQes,eQes,'bs--','LineWidth',1.2,'MarkerSize',5);
        hasQes = true;
    end
end
grid on;
xlabel('Distance (m)');
ylabel('Geometric standard deviation \sigma_g');
title('PSD width evolution');
if hasExp && hasQes
    legend('Experiment (Oak Ridge initial)','QES','Location','best');
elseif hasQes
    legend('QES','Location','best');
elseif hasExp
    legend('Experiment (Oak Ridge initial)','Location','best');
else
    text(0.1,0.5,'No sigma data','Units','normalized');
end
savefig_png(fig,fullfile(outDir,'sigma_evolution.png'));
close(fig);
end

function plot_heatmap(psd,outDir)
distances = unique(psd.qes.distance_m(isfinite(psd.qes.distance_m)));
if isempty(distances)
    return;
end
M = NaN(numel(distances),9);
for k = 1:numel(distances)
    [found,~,distribution] = representative_at_distance(psd.qes,distances(k));
    if found
        M(k,:) = distribution;
    end
end

fig = figure('Color','w');
imagesc(1:numel(distances),1:9,M.');
set(gca,'YTick',1:9,'YTickLabel', ...
    {'1','2','3','4','5','6','7','8','9'});
if numel(distances) <= 20
    set(gca,'XTick',1:numel(distances), ...
        'XTickLabel',cellstr(num2str(distances(:),'%g')));
end
grid on;
xlabel('Distance (m)');
ylabel('PSD bin');
title('QES mass fraction heatmap');
colorbar;
savefig_png(fig,fullfile(outDir,'psd_mass_fraction_heatmap.png'));
close(fig);
end

function [found,distance_m,distribution] = representative_at_distance(qes,target)
found = false;
distance_m = target;
distribution = zeros(1,9);
dist = qes.distance_m(:);
finiteDistance = isfinite(dist);
if ~any(finiteDistance)
    return;
end
uniqueDistance = unique(dist(finiteDistance));
[~,j] = min(abs(uniqueDistance-target));
distance_m = uniqueDistance(j);
rows = find(abs(dist-distance_m) <= max(1.0e-12,1.0e-9*abs(distance_m)));
if isempty(rows)
    return;
end

% The CSV is cumulative within each sampler window, so use the final row
% for each sampler ID before averaging duplicate sampler positions.
ids = qes.sampler_id(rows);
selected = zeros(0,1);
if iscell(ids)
    uid = unique(ids);
    for k = 1:numel(uid)
        q = strcmp(ids,uid{k});
        idx = rows(q);
        [~,jmax] = max(qes.sampling_time_s(idx));
        selected(end+1,1) = idx(jmax); %#ok<AGROW>
    end
else
    [~,jmax] = max(qes.sampling_time_s(rows));
    selected = rows(jmax);
end
distribution = average_rows(qes.mass_fraction(selected,:));
found = any(distribution > 0);
end

function y = average_rows(x)
if isempty(x)
    y = zeros(1,9);
    return;
end
y = zeros(1,size(x,2));
for k = 1:size(x,2)
    v = x(:,k);
    v = v(isfinite(v));
    if isempty(v)
        y(k) = NaN;
    else
        y(k) = mean(v);
    end
end
end

function [distance,value] = grouped_mean(distanceIn,valueIn)
distanceIn = distanceIn(:);
valueIn = valueIn(:);
q = isfinite(distanceIn) & isfinite(valueIn);
distanceIn = distanceIn(q);
valueIn = valueIn(q);
distance = unique(distanceIn);
value = NaN(size(distance));
for k = 1:numel(distance)
    value(k) = mean(valueIn(distanceIn == distance(k)));
end
end
