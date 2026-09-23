function plot_aerosol_psd_results(PSD,metrics,outputDir)
% Plot Oak Ridge/QES initial and QES distance-evolution PSD results.

if nargin < 3 || isempty(outputDir)
    cfg = CESTA1987_config();
    outputDir = fullfile(cfg.outputDir,'PSD_analysis');
end
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end
if ~isfield(PSD,'qes') || isempty(PSD.qes)
    warning('No QES PSD data available for plotting.');
    return;
end

plot_initial_comparison(PSD,outputDir);
plot_distance_evolution(PSD,outputDir);
plot_d50_evolution(metrics,outputDir);
plot_sigma_evolution(metrics,outputDir);
plot_psd_heatmap(PSD,outputDir);
end

function plot_initial_comparison(PSD,outputDir)
qes = PSD.qes;
times = [qes.time];
n = numel(qes);
initialMass = zeros(1,9);
count = 0;
positive = false(1,n);
for i = 1:n
    positive(i) = sum(qes(i).mass_fraction) > 0;
end
if any(positive)
    t0 = min(times(positive));
else
    t0 = min(times);
end
for i = 1:n
    if abs(qes(i).time-t0) <= max(1.0e-12,1.0e-9*abs(t0))
        initialMass = initialMass + qes(i).mass_fraction;
        count = count + 1;
    end
end
if count > 0
    initialMass = initialMass / count;
else
    initialMass(:) = NaN;
end

diameter = representative_diameter_um();
fig = figure('Color','w');
% Oak Ridge input is geometric; the collector bins are aerodynamic. Show
% separate panels instead of a misleading overlay of different quantities.
if isfield(PSD,'exp') && ~isempty(PSD.exp.diameter_um)
    subplot(1,2,1);
    semilogx(PSD.exp.diameter_um,PSD.exp.mass_fraction, ...
        'ko-','LineWidth',1.2,'MarkerSize',5);
    grid on;
    xlabel('Initial geometric diameter (\mum)');
    ylabel('Dry mass fraction');
    title('Oak Ridge initialization reference');
    subplot(1,2,2);
end
semilogx(diameter,initialMass,'bs--','LineWidth',1.2,'MarkerSize',5);
grid on;
xlabel('Aerodynamic observation class (\mum)');
ylabel('Dry UO_2F_2 mass fraction');
title(sprintf('First sampled QES distribution (t=%.3g s)',t0));
savefig_png(fig,fullfile(outputDir,'oakridge_reference_and_qes_sample.png'));
close(fig);
end

function plot_distance_evolution(PSD,outputDir)
targets = [10 20 40 70 100 200 500];
fig = figure('Color','w');
hold on;
labels = {};
plotted = false;
diameter = representative_diameter_um();
for k = 1:numel(targets)
    [found,distance,massFraction] = representative_at_distance(PSD.qes,targets(k));
    if found
        plot(diameter,massFraction,'o-','LineWidth',1.1,'MarkerSize',4);
        labels{end+1} = sprintf('%.0f m',distance); %#ok<AGROW>
        plotted = true;
    end
end
set(gca,'XScale','log');
grid on;
xlabel('Aerodynamic observation class (\mum)');
ylabel('Dry UO_2F_2 mass fraction');
title('QES PSD evolution by distance');
if plotted
    legend(labels,'Location','best');
else
    text(0.1,0.5,'No QES PSD rows','Units','normalized');
    axis off;
end
savefig_png(fig,fullfile(outputDir,'PSD_evolution_distance.png'));
close(fig);
end

function plot_d50_evolution(metrics,outputDir)
fig = figure('Color','w');
hold on;
[distance,d50] = grouped_mean(metrics.distance_m,metrics.d50_um);
if any(isfinite(d50))
    plot(distance,d50,'bs-','LineWidth',1.2,'MarkerSize',5);
end
grid on;
xlabel('Distance (m)');
ylabel('Mass median aerodynamic diameter (\mum)');
title('QES d50 evolution');
if any(isfinite(d50))
    legend('QES','Location','best');
else
    text(0.1,0.5,'No d50 data','Units','normalized');
end
savefig_png(fig,fullfile(outputDir,'d50_evolution.png'));
close(fig);
end

function plot_sigma_evolution(metrics,outputDir)
fig = figure('Color','w');
hold on;
[distance,sigma] = grouped_mean(metrics.distance_m,metrics.sigma_g);
if any(isfinite(sigma))
    plot(distance,sigma,'bs-','LineWidth',1.2,'MarkerSize',5);
end
grid on;
xlabel('Distance (m)');
ylabel('Geometric standard deviation \sigma_g');
title('QES PSD width evolution');
if any(isfinite(sigma))
    legend('QES','Location','best');
else
    text(0.1,0.5,'No sigma data','Units','normalized');
end
savefig_png(fig,fullfile(outputDir,'sigma_evolution.png'));
close(fig);
end

function plot_psd_heatmap(PSD,outputDir)
distances = unique([PSD.qes.distance]);
distances = distances(isfinite(distances));
if isempty(distances)
    return;
end
M = NaN(numel(distances),9);
for k = 1:numel(distances)
    [found,~,massFraction] = representative_at_distance(PSD.qes,distances(k));
    if found
        M(k,:) = massFraction;
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
ylabel('Aerodynamic observation bin');
title('QES mass fraction');
colorbar;
savefig_png(fig,fullfile(outputDir,'PSD_heatmap.png'));
close(fig);
end

function [found,distance,massFraction] = representative_at_distance(qes,target)
found = false;
distance = target;
massFraction = zeros(1,9);
distances = [qes.distance];
finiteDistance = isfinite(distances);
if ~any(finiteDistance)
    return;
end
uniqueDistance = unique(distances(finiteDistance));
[~,j] = min(abs(uniqueDistance-target));
distance = uniqueDistance(j);
rows = find(abs(distances-distance) <= max(1.0e-12,1.0e-9*abs(distance)));
if isempty(rows)
    return;
end

ids = {qes(rows).sampler_id};
selected = zeros(0,1);
if iscell(ids) && ~isempty(ids)
    uid = unique(ids);
    for k = 1:numel(uid)
        group = rows(strcmp(ids,uid{k}));
        times = [qes(group).time];
        [~,jmax] = max(times);
        selected(end+1,1) = group(jmax); %#ok<AGROW>
    end
else
    times = [qes(rows).time];
    [~,jmax] = max(times);
    selected = rows(jmax);
end

mass = zeros(numel(selected),9);
for k = 1:numel(selected)
    mass(k,:) = qes(selected(k)).mass_fraction;
end
massFraction = mean_finite_rows(mass);
found = any(massFraction > 0);
end

function y = mean_finite_rows(x)
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

function diameter = representative_diameter_um()
diameter = [0.2 sqrt(0.4*0.7) sqrt(0.7*1.1) sqrt(1.1*2.1) ...
    sqrt(2.1*3.3) sqrt(3.3*4.7) sqrt(4.7*5.8) sqrt(5.8*9.0) 12.0];
end
